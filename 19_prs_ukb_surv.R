## Survival analysis

library(tidyverse)
library(survival)

reg_data <- rio::import("../files/reg_data.RData")

prs_in_bmi <- rio::import("../files/prs_in_bmi.tsv")

causes <- paste(c("all", "cv", "t2d"), "cause", sep = "_")

cox_res <- lapply(
    causes,
    function(cause){
        cox_formula <- formula(paste(paste0("Surv(time, ", cause, ")"),
                                     paste("prs_conc + prs_disc", reg_data$covariate_term,
                                           sep = " + "),
                                     sep  = " ~ "))
        return(coxph(cox_formula, data = reg_data$reg_df))
    }
)

cox_tidy <- lapply(cox_res,
                   function(x) {
                       s <- summary(x)
                       cbind(s$coefficients[1:2,],
                             s$conf.int[1:2,3:4]) %>%
                           data.frame %>%
                           rownames_to_column %>%
                           `names<-`(c("profile", "b", "hr", "se",
                                       "z", "p", "ci_lo", "ci_up"))
                   }) %>%
    `names<-`(causes) %>%
    bind_rows(.id = "cause") %>%
    inner_join(prs_in_bmi, by = c("profile" = "term")) %>%
    mutate(across(c(b, se), ~ .x / u)) %>%
    transmute(cause, profile,
              hr = exp(b),
              ci_lo = exp(b - qnorm(1-0.05/2)*se),
              ci_up = exp(b + qnorm(1-0.05/2)*se),
              p, units = "hr per 1kg/m2")

rio::export(cox_tidy, "../files/cox_tidy.tsv")
