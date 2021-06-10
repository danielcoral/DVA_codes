## MR Plot

library(tidyverse)
library(ggh4x)

## Results from all analysis

mr_res <- rio::import("~/dva/files/mr_res.tsv")

trait_order <- c("Urate", "WHR", "#Leuco",
                 "SBP", "DBP",
                 "HDL", "TG", "LDL",
                 "FVC",
                 "HeelBMD", "Phosphate",
                 "TotalLeanMass",
                 "#Retic", "Platelet volume",
                 "ArmFatMassR", "LegFatMassR")

trait_clusters <- tribble(
    ~ trait_short, ~ clus,
    "Urate", 1,
    "WHR", 1,
    "#Leuco", 1,
    "SBP", 2,
    "DBP", 2,
    "HDL", 3,
    "TG", 3,
    "LDL", 3,
    "FVC", 4,
    "HeelBMD", 5,
    "Phosphate", 5,
    "TotalLeanMass", 6,
    "#Retic", 7,
    "Platelet volume", 7,
    "ArmFatMassR", 8,
    "LegFatMassR", 8
)

sig_mr <- mr_res %>%
    filter(method == "mr_ivw") %>%
    mutate(padj = p.adjust(pval, "fdr")) %>%
    filter(padj < 0.05)

panelsizes <- trait_clusters %>%
    group_by(clus) %>%
    summarise(n = n())

mr_res %>%
    filter(!(method == "conmix" & psi != 1.5)) %>%
    inner_join(trait_clusters) %>%
    mutate(exposure = factor(trait_short, levels = rev(trait_order)),
           sig = ifelse("|"(method == "mr_ivw" & trait_short %in% sig_mr$trait_short,
                            method == "intercept" & pval < 0.05),
                        23, 1),
           method = recode_factor(method,
                                  "mr_ivw" = "IVW \u03b2", "intercept" = "Egger \u03b2\u2080",
                                  "egger" = "Egger \u03b2", "wmedian" = "W median \u03b2",
                                  "conmix" = "ConMix \u03b2")) %>%
    ggplot(aes(b, exposure, shape = sig)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_errorbar(aes(xmin = ci_lo, xmax = ci_up),
                  alpha = .4, size = .4, width = .2) +
    geom_point(aes(fill = ifelse(sig == 23, "black", "grey"))) +
    scale_shape_identity() + scale_fill_identity() +
    facet_grid(clus ~ method, scales = "free") +
    force_panelsizes(rows = panelsizes$n) +
    labs(x = "LogOR per SD") +
    theme_light() +
    theme(legend.position = "none",
          axis.title.y = element_blank(),
          strip.background.y = element_blank(),
          strip.text.y = element_blank(),
          axis.text.x = element_text(size = 5.5))

ggsave("../plots/mr_res_all.png", width = 170, units = "mm")

## ConMix across different psi
mr_res %>%
    filter(method == "conmix") %>%
    inner_join(trait_clusters) %>%
    mutate(exposure = factor(trait_short, levels = rev(trait_order))) %>%
    ggplot(aes("", b)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_up),
                  alpha = .4, size = .4, width = .2) +
    geom_point() +
    facet_grid(exposure ~ psi, scales = "free", switch = "both") +
    labs(y = "LogOR per SD\n", x = "\u03a8 multiplier") +
    scale_y_continuous(position = "right") +
    theme_light() +
    theme(strip.background = element_blank(),
          strip.text = element_text(color = "black"),
          strip.text.y.left = element_text(angle = 0),
          axis.ticks.y = element_blank(),
          axis.text.y = element_text(size = 5))

ggsave("../plots/mr_res_conmix.png", height = 225, width = 170, units = "mm")
