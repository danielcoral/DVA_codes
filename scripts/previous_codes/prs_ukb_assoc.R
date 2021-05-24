#### Scores associations in UK Biobank

library(tidyverse)

reg_data <- rio::import("../files/reg_data.tsv")

## Results will be in SD units per 1SD of PRS
reg_lm <- map(select(reg_data, -c(eid:t2d_dx0, conc_prs:time)),
               ~ lm(.x ~ conc_prs + disc_prs, data = reg_data) %>%
                   broom::tidy(conf.int = TRUE) %>%
                   slice(2:3)) %>%
    bind_rows(.id = "trait") %>%
    mutate(units = "sd/prs_sd")

reg_t2d <- glm(t2d_dx0 ~ conc_prs + disc_prs +
                   age0 + sex +
                   genpcs_1 + genpcs_2 + genpcs_3 + genpcs_4 + genpcs_5 +
                   genpcs_6 + genpcs_7 + genpcs_8 + genpcs_9 + genpcs_10,
               data = reg_data, family = "binomial") %>%
    {
        bind_cols(broom::tidy(.)[2:3,],
                  confint.default(.)[2:3,] %>%
                  data.frame %>%
                  `names<-`(c("conf.low", "conf.high")))
    } %>%
    mutate(trait = "t2d",
           units = "logOR/prs_sd")

reg_res <- bind_rows(reg_lm, reg_t2d)

## SD of BMI in UKB
bmi_sd <- rio::import("../files/ukb_f.RData") %>%
    pull(bmi0) %>%
    sd(na.rm = T)

## 1 SD of each score in BMI units
prs_in_bmi <- reg_res %>%
    filter(trait == "bmi0") %>%
    transmute(term, u = estimate * bmi_sd)

rio::export(prs_in_bmi, "../files/prs_in_bmi.tsv")

## Expressing betas in terms of change per BMI units predicted by each PRS
ukb_prs_res <- reg_res %>%
    inner_join(prs_in_bmi, by = "term") %>%
    mutate(units = ifelse(trait == "bmi0", units,
                          gsub("/prs_sd", "/1kg_m2", units)),
           across(c(estimate, std.error, conf.low, conf.high),
                  ~ ifelse(trait != "bmi0", .x / u, .x)),
           u = NULL) %>%
    ## Comparing effects
    group_by(trait) %>%
    mutate(bcomp = abs(estimate[1] - estimate[2]),
           secomp = sqrt(std.error[1]^2 + std.error[2]^2),
           pcomp = pnorm(-bcomp/secomp)) %>%
    ungroup %>%
    mutate(padj = p.adjust(pcomp, "fdr"))

rio::export(ukb_prs_res, "../files/ukb_prs_res.tsv")
