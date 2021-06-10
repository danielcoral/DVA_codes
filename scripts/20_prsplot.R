## Combining PRS results

library(tidyverse)
library(patchwork)
library(ggh4x)

p <- rio::import("~/dva/files/pwbiovu_plot.RData")

l <- rio::import("~/dva/files/lwbiovu_plot.RData") +
    theme(axis.title.x = element_text(size = 7), plot.margin = unit(c(0,.05, 0.1, .05), "cm"))

s <- rio::import("~/dva/files/prs_surv.RData") +
    theme(axis.title.x = element_text(size = 7),
          axis.text.y = element_text(size = 7),
          plot.margin = unit(c(0,.05, 0.1, .05), "cm"))

wrap_plots(p, wrap_plots(l, s, ncol = 1, heights = c(.6, .4), tag_level = "keep")) +
    plot_annotation(tag_levels = "A")

ggsave("../plots/prsplot.png", width = 170, units = "mm")
