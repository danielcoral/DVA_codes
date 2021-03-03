## Plotting profiles

library(tidyverse)

mix <- rio::import("../files/mix.tsv")

prof <- mix %>%
    mutate(beta.t2d = exp(beta.t2d)) %>%
    ggplot(aes(beta.bmi, beta.t2d, colour = factor(disc))) +
    geom_point(size = 3, alpha = 0.5) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_colour_manual(values = c("red", "blue")) +
    geom_label(aes(max(beta.bmi), 1, label = "Concordant SNPs"),
               label.padding = unit(0.1, "cm"),
               hjust = 1, vjust = -1.01, colour = "red", size = 6) +
    geom_label(aes(max(beta.bmi), 1, label = "Discordant SNPs"),
               label.padding = unit(0.1, "cm"),
               hjust = 1, vjust = 2.01, colour = "blue", size = 6) +
    theme(legend.position = "none",
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)) +
    labs(x = "BMI increase (SD units)", y = "T2D OR")

ggsave("../plots/profiles.png", prof)
