## MR Plot

library(tidyverse)

mr_res_univariate <- rio::import("~/dva/files/mr_res_univariate.tsv")

mr_profiles <- rio::import("~/dva/files/mr_profiles.tsv")

trait_order <- c("SHBG", "Urate", "GGT", "AST", "ALT", "WHR", "AdipoQ", "#Leuco",
                 "SBP", "DBP",
                 "HDL", "TG", "LDL", "ApoA1",
                 "FVC",
                 "HeelBMD", "Phosphate",
                 "TotalLeanMass",
                 "#Retic",
                 "ArmFatMassR", "LegFatMassR")

mr_uni <- mr_res_univariate %>%
    filter(method == "ivw") %>%
    transmute(b, ci_lo, ci_up,
              exposure = factor(exposure, levels = rev(trait_order)),
              analysis = recode(analysis, all = "All instruments",
                                bmi_sig = "BMI significant"))

mr_p <- mr_profiles %>%
    filter(method == "ivw") %>%
    transmute(b, ci_lo, ci_up,
              exposure = factor(exposure, levels = rev(trait_order)),
              analysis = disc)

mr_res <- bind_rows(mr_uni, mr_p) %>%
    ggplot(aes(b, exposure)) +
    geom_point() +
    geom_errorbar(aes(xmin = ci_lo, xmax = ci_up),
                  alpha = .4, size = .4, width = .2) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    facet_wrap(~ analysis, nrow = 1) +
    labs(x = "LogOR per SD", y = NULL) +
    theme_light()

ggsave("../plots/mr_res.png")
