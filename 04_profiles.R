## Plotting profiles

library(tidyverse)

mix <- rio::import("../files/mix.txt")

prof <- ggplot(mix, aes(beta.bmi, beta.t2d, colour = factor(disc))) +
    geom_point(size = 3, alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_colour_manual(values = c("red", "blue")) +
    geom_label(aes(max(beta.bmi), 0, label = "Concordant"),
               label.padding = unit(0.1, "cm"),
               hjust = 1, vjust = -1.01, colour = "red", size = 6) +
    geom_label(aes(max(beta.bmi), 0, label = "Discordant"),
               label.padding = unit(0.1, "cm"),
               hjust = 1, vjust = 2.01, colour = "blue", size = 6) +
    theme(legend.position = "none",
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)) +
    labs(x = "Beta BMI", y = "Beta T2D")

ggsave("../plots/profiles.png", prof)
