## Merging phenotypes with PRS and mortality - UKB

library(tidyverse)

## Phenotype data
ukb_f <- rio::import("../files/ukb_f.RData")

## Function for INT (using Blom's offset)
INT <- function(x, k = 3/8){
    qnorm((rank(x, na.last = "keep") - k) / ((sum(!is.na(x))) - (2 * k) + 1))
}

## Residualization and INT of continuous variables
ukb_f <- ukb_f %>%
    mutate(across(c(bmi0:ldl0),
                  ~INT(resid(lm(.x ~ age0 + sex +
                                    genpcs_1 + genpcs_2 + genpcs_3 + genpcs_4 + genpcs_5 +
                                    genpcs_6 + genpcs_7 + genpcs_8 + genpcs_9 + genpcs_10,
                                na.action = na.exclude)))
                  ))

## Genetic scores
all_scores <- list(
    rio::import("../files/conc_prs.tsv") %>% `names<-`(c("eid", "conc_prs")),
    rio::import("../files/disc_prs.tsv") %>% `names<-`(c("eid", "disc_prs"))
) %>%
    reduce(inner_join) %>%
    mutate(across(-eid, scale))

## Mortality data (Extracted manually from UKB data portal)
mor <- rio::import("../files/death.txt")
morc <- rio::import("../files/death_cause.txt")
mor_ed <- mor %>%
    filter(ins_index == 0) %>%
    transmute(eid,
              date_of_death = as.Date(format(as.Date(date_of_death,
                                                     format = "%d/%m/%Y"),
                                             "%Y-%m-%d")),
              cv_m = eid %in% unique(morc$eid[grep("I", morc$cause_icd10)]),
              t2d_m = eid %in% unique(morc$eid[grep("E11", morc$cause_icd10)]))

## Merging data 
reg_data <- inner_join(ukb_f, all_scores, by = "eid") %>% ## Phenotype data with scores
    left_join(mor_ed, by = "eid") %>%                   ## ...and then with mortality data
    ## Parsing mortality information
    mutate(death_bf_limit = date_of_death < "2020-07-31",
           death_recorded = !is.na(date_of_death),
           all_cause = ifelse(death_recorded & death_bf_limit, 1, 0),
           cv_cause = ifelse(all_cause == 1 & cv_m, 1, 0),
           t2d_cause = ifelse(death_recorded & death_bf_limit & t2d_m, 1, 0),
           time = round(ifelse(all_cause == 1,
                               as.numeric((date_of_death - date0)),
                               as.numeric((as.Date("2020-07-31") - date0))) / 30.44))

rio::export(reg_data, "../files/reg_data.tsv")
