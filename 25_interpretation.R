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

## Percent change in odds of diseases
ukb_prs %>%
    filter(grepl("logOR", units)) %>%
    transmute(trait, term,
              across(c(estimate, lo, up),
                     ~ 100 * (exp(.x) - 1)),
              mantissa, exponent)

## Change of continuous traits in natural units
std_devs <- rio::import("../files/ukb_f.RData") %>%
    summarise(across(all_of(ukb_prs %>%
                            filter(grepl("^sd", units)) %>%
                            pull(trait)),
                     sd, na.rm = T)) %>%
    pivot_longer(everything(), names_to = "trait", values_to = "sd_value")

ukb_prs %>%
    filter(units == "sd/1kg_m2") %>%
    left_join(std_devs, by = "trait") %>%
    transmute(trait, term,
              eff = estimate * sd_value,
              eff = ifelse(grepl("fat", trait), eff * 1000, eff),
              units = case_when(
                  trait == "whr" ~ "ratio unit",
                  grepl("fat", trait) ~ "grams",
                  trait %in% c("sbp","dbp") ~ "mmHg",
                  trait %in% c("hdl", "ldl", "tg", "t_chol", "ldl") ~ "mmol/L"
              ),
              units = paste0(units, " / 1kg/m2"),
              effround = round(eff, 2),
              mantissa, exponent) %>%
    group_split(term) %>%
    lapply(data.frame)

## Change of SD per SD of PRS in fat mass
prs_in_bmi <- rio::import("../files/prs_in_bmi.tsv")
ukb_prs %>%
    filter(grepl("fat", trait)) %>%
    left_join(prs_in_bmi, by = "term") %>%
    left_join(std_devs, by = "trait") %>%
    transmute(trait, term,
              estimate = estimate * sd_value * u * 1000,
              units = "grams / 1SD PRS", mantissa, exponent) %>%
    group_split(term) %>%
    lapply(data.frame)
    

## Survival analysis
cox_tidy <- rio::import("../files/cox_tidy.tsv")

cox_tidy

## BioVU
pw <- rio::import("../files/pw_biovu.tsv")
lw <- rio::import("../files/lw_biovu.tsv")

pw %>% filter(description == "Obesity") %>%
    mutate(across(c(beta, ci_up, ci_lo), exp))

pw %>% filter(description == "Type 2 diabetes") %>%
    mutate(across(c(beta, ci_up, ci_lo), exp))

pw %>% filter(description == "Type 1 diabetes") %>%
    mutate(across(c(beta, ci_up, ci_lo), exp))

pw %>% filter(description == "Essential hypertension") %>%
    mutate(across(c(beta, ci_up, ci_lo), exp))
