## Survival analysis with PRS

library(tidyverse)
library(survival)
library(ggh4x)

ukbphen <- rio::import("~/dva/files/ukbphen.RData")
    
## Formulas for Cox regression
causes <- paste(c("all", "cv", "t2d"), "cause", sep = "_")
covars <- paste("age0", "sex",
                paste0("gpc", 1:10, collapse = " + "),
                sep = " + ")
predictors <- c("conc", "disc",
                str_subset(names(ukbphen), "snpclus_"))

## Running the formulas
cox_res <- map(causes,
               function(cause){
                   survterm <- paste0("Surv(time, ", cause, ")")
                   map(predictors,
                       function(p){
                           f <- paste0(survterm, " ~ ", p, " + ", covars)
                           coxph(formula(f), data = ukbphen)
                       }) %>%
                       `names<-`(predictors)
               }) %>%
    `names<-`(causes)

## Tidying up results
cox_tidy <- map(cox_res,
                function(cause_res){
                    map(cause_res,
                        function(mod){
                            summary(mod) %>%
                                "$"(coefficients) %>%
                                data.frame %>%
                                slice(1) %>%
                                rownames_to_column %>%
                                `names<-`(c("predictor", "b", "hr", "se", "z", "p"))
                        }) %>%
                        bind_rows()
                }) %>%
    bind_rows(.id = "cause") %>%
    mutate(ci_lo = exp(b - qnorm(1-0.05/2)*se),
           ci_up = exp(b + qnorm(1-0.05/2)*se))

rio::export(cox_tidy, "~/dva/files/cox_tidy.tsv")

## Plotting concordant and discordant effects
mix <- rio::import("~/dva/files/mix.tsv")

nsnps_prs <- mix %>%
    mutate(predictor = ifelse(disc == 1, "disc", "conc")) %>%
    group_by(predictor) %>%
    summarise(nsnps = n())

survplot_dat <- cox_tidy %>%
    mutate(cause = recode_factor(cause,
                                 all_cause = "All cause",
                                 cv_cause = "Cardiovascular",
                                 t2d_cause = "Diabetes complications"),
           predictor_groups = factor(case_when(
               predictor %in% c("conc", "disc") ~ "Profiles",
               grepl("snpclus_conc", predictor) ~ "Concordant\nsubgroups",
               grepl("snpclus_disc", predictor) ~ "Discordant\nsubgroups"
           ), levels = c("Profiles", paste(c("Concordant", "Discordant"),
                                           "subgroups", sep = "\n")))) %>%
    left_join(nsnps_prs) %>%
    mutate(nsnps = coalesce(nsnps, as.numeric(gsub("snpclus_.+_n", "", predictor))),
           predictor = str_replace_all(predictor,
                                       c("snpclus_.+_n" = "N SNPs = ",
                                         "conc" = "Concordant",
                                         "disc" = "Discordant")))

prs_surv <- survplot_dat %>%
    ggplot(aes(hr, reorder(predictor, nsnps))) +
    geom_errorbar(aes(xmin = ci_lo, xmax = ci_up),
                  alpha = .4, size = .4, width = .2) +
    geom_point(aes(color = predictor == "Discordant" | grepl("Discordant", predictor_groups))) +
    geom_point(shape = 1, col = "black", alpha = .4) +
    scale_colour_discrete(labels = c("Concordant", "Discordant")) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    labs(title = "", x = "HR per allele") +
    facet_grid(predictor_groups ~ cause, scales = "free_y", switch = "y") +
    theme_light() +
    theme(axis.title.y = element_blank(),
          strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, face = "bold"),
          legend.title = element_blank(), legend.position = "top") +
    force_panelsizes(rows = c(.5, .25, .25))

ggsave("../docs/plots/prs_surv.png")
