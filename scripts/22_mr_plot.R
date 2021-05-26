## MR Plot

library(tidyverse)
library(ggh4x)

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

mr_res_dat <- bind_rows(mr_uni, mr_p) %>%
    mutate(clus = case_when(
               exposure %in% c("SHBG", "Urate", "GGT", "AST", "ALT",
                               "WHR", "AdipoQ", "#Leuco") ~ 1,
               exposure %in% c("SBP", "DBP") ~ 2,
               exposure %in% c("HDL", "TG", "LDL", "ApoA1") ~ 3,
               exposure == "FVC" ~ 4,
               exposure %in% c("HeelBMD", "Phosphate") ~ 5,
               exposure == "TotalLeanMass" ~ 6,
               exposure == "#Retic" ~ 7,
               exposure %in% c("ArmFatMassR", "LegFatMassR") ~ 8
           ),
           strat = ifelse(analysis %in% c("Concordant", "Discordant"),
                          "Stratified by T2D Beta", analysis),
           analysis = ifelse(analysis %in% c("All instruments", "BMI significant"),
                             "", analysis))

panelsizes <- mr_res_dat %>%
    group_by(clus) %>%
    summarise(n = length(unique(exposure)))

mr_res_plot <- mr_res_dat %>%
    ggplot(aes(b, exposure)) +
    geom_point() +
    geom_errorbar(aes(xmin = ci_lo, xmax = ci_up),
                  alpha = .4, size = .4, width = .2) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    facet_nested(clus ~ strat + analysis, scales = "free_y") +
    labs(x = "LogOR per SD", y = NULL) +
    theme_light() +
    theme(strip.background.y = element_blank(),
          strip.text.y = element_blank()) +
    force_panelsizes(rows = panelsizes$n,
                     cols = c(.3,.3,.2,.2))

ggsave("../docs/plots/mr_res.png", width = 10)
