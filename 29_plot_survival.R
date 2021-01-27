## Plot of effects of concordant and discordant profiles on survival
library(tidyverse)

cox_tidy <- rio::import("../files/cox_tidy.tsv")

survival_plot <- cox_tidy %>%
    mutate(cause = recode_factor(cause,
                                 all_cause = "All",
                                 cv_cause = "CV",
                                 t2d_cause = "Diabetes")) %>%
    ggplot(aes(hr, profile)) +
    geom_errorbar(aes(xmin = ci_lo, xmax = ci_up),
                  alpha = .4, size = .4, width = .2) +
    geom_point(aes(colour = profile), alpha = .5) +
    geom_point(shape = 1, col = "black", alpha = .4) +
    scale_colour_manual(values = c("red", "blue"),
                        labels = c("Concordant", "Discordant"),
                        guide = guide_legend(title = NULL, direction = "vertical")) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    labs(title = "", x = "HR per 1kg/m\u00B2") +
    facet_wrap(cause~., strip.position = "top", ncol = 1) +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(size = 8),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "bottom",
          strip.text.y.left = element_text(angle = 0, margin = margin()))
