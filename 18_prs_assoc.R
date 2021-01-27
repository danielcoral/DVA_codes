#### Scores associations in UK Biobank

library(tidyverse)

## Data for regression
reg_data <- rio::import("../files/reg_data.RData")

## Applying the formulas
reg_res <- lapply(
    seq_along(reg_data$endpoints),
    function(i){
        outcome <- reg_data$endpoints[i]
        f <- paste(paste(outcome, reg_data$scores_term, sep = " ~ "),
                   reg_data$covariate_term, sep = " + ")
        fam <- ifelse(reg_data$endpoint_types[i] == "numeric", "gaussian", "binomial")
        glm(formula = as.formula(f), family = fam, data = reg_data$reg_df)
    }
)

ukb_prs_res <- lapply(
    reg_res,
    function(x)
        bind_cols(broom::tidy(x)[2:3,],
                  setNames(data.frame(confint.default(x)[2:3,],
                                      row.names = NULL),
                           c("lo", "up")),
                  units = x$family$family)
) %>%
    `names<-`(reg_data$endpoints) %>%
    bind_rows(.id = "trait") %>%
    ## Units of effects
    mutate(units = ifelse(units == "gaussian", "sd/prs_sd", "logOR/prs_sd"))

## SD of BMI in UKB
bmi_sd <- rio::import("../files/ukb_f.RData") %>%
    pull(bmi) %>%
    sd(na.rm = T)

## 1 SD of each score in BMI units
prs_in_bmi <- ukb_prs_res %>%
    filter(trait == "bmi") %>%
    transmute(term, u = estimate * bmi_sd)

rio::export(prs_in_bmi, "../files/prs_in_bmi.tsv")

## Expressing betas in terms of change per BMI units predicted by each PRS
ukb_prs_res <- ukb_prs_res %>%
    inner_join(prs_in_bmi, by = "term") %>%
    mutate(units = ifelse(trait == "bmi", units,
                          gsub("/prs_sd", "/1kg_m2", units)),
           across(c(estimate, std.error, lo, up),
                  ~ ifelse(trait != "bmi", .x / u, .x)),
           u = NULL)

rio::export(ukb_prs_res, "../files/ukb_prs_res.tsv")
