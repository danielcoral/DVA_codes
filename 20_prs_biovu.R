## Results from BioVU

library(tidyverse)

prs_in_bmi <- rio::import("../files/prs_in_bmi.tsv")

## PheWAS
pw <- vroom::vroom(list.files("../files/biovu/phewas", pattern = "w_unit", full.names = T),
                   id = "profile") %>%
    drop_na %>%
    mutate(profile = gsub(".*phewas_(.*)ordant_GRS.*", "prs_\\1", profile)) %>%
    left_join(prs_in_bmi, by = c("profile" = "term")) %>%
    mutate(prs_sd = log(OR.per.sd) / beta,
           beta = beta * prs_sd / u, SE = SE * prs_sd / u,
           ci_up = beta + qnorm(1-0.05/2)*SE, ci_lo = beta - qnorm(1-0.05/2)*SE,
           priority = description %in% c("Obesity", "Type 2 diabetes",
                                         "Essential hypertension", "Hyperlipidemia",
                                         "Coronary atherosclerosis")) %>%
    group_by(profile) %>%
    mutate(Bonferroni = p.adjust(p, "bonferroni") < 0.05) %>%
    { ## Significant signals
        df <- ungroup(.)
        sig_dx <- df %>%
            filter(Bonferroni) %>%
            pull(phecode) %>%
            unique
        filter(df, priority | phecode %in% sig_dx) } %>%
    transmute(dataset = "all", profile, phecode, description,
              beta, SE, ci_up, ci_lo, p) %>%
    { ## Excluding individuals with T1D
        df <- .
        not1d <- vroom::vroom(list.files("../files/biovu/phewas_rm_t1d",
                                pattern = "w_unit", full.names = T),
                     id = "profile") %>%
            filter(phecode %in% df$phecode) %>%
            mutate(profile = gsub(".*phewas_(.*)ordant_GRS.*", "prs_\\1", profile)) %>%
            left_join(prs_in_bmi, by = c("profile" = "term")) %>%
            mutate(prs_sd = log(OR.per.sd) / beta,
                   beta = beta * prs_sd / u, SE = SE * prs_sd / u,
                   ci_up = beta + qnorm(1-0.05/2)*SE, ci_lo = beta - qnorm(1-0.05/2)*SE) %>%
            transmute(dataset = "no_t1d", profile, phecode, description,
                      beta, SE, ci_up, ci_lo, p)
        bind_rows(df, not1d) }

rio::export(pw, "../files/pw_biovu.tsv")

## LabWAS
lw <- vroom::vroom(list.files("../files/biovu/labwas", full.names = T), id = "path") %>%
    filter(Predictor == "Unit_PRS") %>% ## Unweighted PRS
    mutate(pop_profile = gsub(".*labwas/(.*)_PRS.*", "\\1", path), path = NULL) %>%
    separate(pop_profile, c("pop", "profile")) %>%
    mutate(profile = tolower(paste0("prs_", gsub("ordant", "", profile)))) %>%
    left_join(prs_in_bmi, by = c("profile" = "term")) %>%
    mutate(beta = Estimate / u, SE = SE / u,
           ci_up = beta + qnorm(1-0.05/2)*SE, ci_lo = beta - qnorm(1-0.05/2)*SE) %>%
    { ## Significant signals
        df <- .
        sig_lab_EUR <- df %>%
            filter(pop == "EA" & Bonferroni) %>%
            pull(Lab) %>%
            unique
        sig_lab_AFR <- df %>%
            filter(pop == "AA" & Bonferroni) %>%
            pull(Lab) %>%
            unique
        df %>%
            filter(pop == "EA" & Lab %in% sig_lab_EUR | pop == "AA" & Lab %in% sig_lab_AFR) } %>%
    select(pop, profile, Full_name, beta, SE, ci_up, ci_lo, p)

rio::export(lw, "../files/lw_biovu.tsv")
