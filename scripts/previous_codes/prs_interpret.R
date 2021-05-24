## Interpretation of results

library(tidyverse)

## PRS results
ukb_prs <- rio::import("../files/ukb_prs_res.tsv")

## Formatting p-values
ukb_prs <- ukb_prs %>%
    mutate(p = (log(2) +
                pnorm(abs(statistic), lower.tail = F, log.p = T)) /
               log(10),
           mantissa = round(10^(p %% 1), 2),
           exponent = p %/% 1,
           p = NULL)

## T2D OR
ukb_prs %>%
    filter(grepl("logOR", units)) %>%
    transmute(trait, term,
              across(c(estimate, conf.low, conf.high), exp),
              mantissa, exponent)

## Change of continuous traits in natural units
traits <- ukb_prs %>%
    filter(units == "sd/1kg_m2") %>%
    pull(trait) %>%
    unique

traits

nat_u <- c("ratio unit", rep("grams", 2),
           rep("mmHg", 2), rep("mmol/L", 3)) %>%
    `names<-`(traits)

std_devs <- rio::import("../files/ukb_f.RData") %>%
    summarise(across(all_of(traits), sd, na.rm = T)) %>%
    pivot_longer(everything(), names_to = "trait", values_to = "sd_value") %>%
    mutate(sd_value = ifelse(trait %in% c("armfat0", "legfat0"),
                             sd_value * 1000, sd_value),
           nat.units = str_replace_all(trait, nat_u))

ukb_prs %>%
    right_join(std_devs, by = "trait") %>%
    transmute(trait, term,
              across(c(estimate, conf.low, conf.high),
                     ~ round(.x * sd_value, 3)),
              units = paste0(nat.units, " / 1kg/m2"),
              mantissa, exponent) %>%
    group_split(term) %>%
    lapply(data.frame)

## Change of SD per SD of PRS in adiposity traits
prs_in_bmi <- rio::import("../files/prs_in_bmi.tsv")

ukb_prs %>%
    left_join(prs_in_bmi, by = "term") %>%
    right_join(std_devs, by = "trait") %>%
    filter(grepl("fat|vol", trait)) %>%
    transmute(trait, term,
              across(c(estimate, conf.low, conf.high),
                     ~ (sd.u = round(.x * u, 3))),
              units = paste0("SD / PRS SD"),
              mantissa, exponent) %>%    
    group_split(term) %>%
    lapply(data.frame)

## Change of natural units per SD of PRS in adiposity traits
ukb_prs %>%
    left_join(prs_in_bmi, by = "term") %>%
    right_join(std_devs, by = "trait") %>%
    filter(grepl("fat|vol", trait)) %>%
    transmute(trait, term,
              across(c(estimate, conf.low, conf.high),
                     ~ round(.x * sd_value * u, 3)),
              units = paste0(nat.units, " / PRS SD"),
              mantissa, exponent) %>%    
    group_split(term) %>%
    lapply(data.frame)

## Survival analysis
cox_tidy <- rio::import("../files/cox_tidy.tsv")

cox_tidy

## BioVU
pw <- rio::import("../files/pw_biovu.tsv")
lw <- rio::import("../files/lw_biovu.tsv")

pw %>% filter(description == "Obesity") %>%
    mutate(across(c(beta, conf.high, conf.low), exp))

pw %>% filter(description == "Type 2 diabetes") %>%
    mutate(across(c(beta, conf.high, conf.low), exp))

pw %>% filter(description == "Type 1 diabetes") %>%
    mutate(across(c(beta, conf.high, conf.low), exp))

pw %>% filter(description == "Essential hypertension") %>%
    mutate(across(c(beta, conf.high, conf.low), exp))
