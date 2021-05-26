## Analysis pipeline plot

library(tidyverse)
library(ggrepel)

aline_dat <- data.frame(
    items = c(NA,
              "Cross-reference\nBMI and T2D\nGWAS",
              "Concordant\nand\ndiscordant\nobesity\nsubtypes",
              "Phenome-wide\nscan",
              "Phenome-wide\ncomparison",
              "Feature\nselection",
              "External\nvalidation\nin BioVU",
              "Subtype-specific\nmortality in UKB",
              "Mendelian\nrandomization"),
    x = seq(0, 8, 1), y = 0,
    ylabels = c(0, rep(c(1,-1), 4))
)

aline_plot <- aline_dat %>% 
    ggplot(aes(x, y)) +
    geom_segment(aes(x = 0, y = 0, xend = max(x) + 1, yend = 0),
                 arrow = arrow(length = unit(0.5, "cm"))) +
    geom_segment(aes(x = x, y = 0, xend = x, yend = ylabels)) +
    geom_point(aes(fill = items), size = 5,
               shape = 21, color = "black") +
    geom_label(aes(y = ylabels, label = items)) +
    guides(fill = "none") +
    ylim(-1.5, 1.5) +
    theme_void()

save(aline_plot, file = "../plot_files/aline_plot.RData")

ggsave("../docs/plots/aline_plot.png", aline_plot)
