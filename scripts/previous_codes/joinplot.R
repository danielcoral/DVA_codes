## Joining MR and PRS plots

lapply(list.files(pattern = "prsplot"), source)

library(cowplot)

joinplot <- plot_grid(plot_grid(prs_plot, surv_plot, ncol = 2,
                                rel_widths = c(.6, .4), labels = c("A", "B")),
                      plot_grid(pwbiovu_plot, lwbiovu_plot,
                                ncol = 2, rel_widths = c(.65,.35), labels = c("C", "D")),
                      ncol = 1, rel_heights = c(.25,.75))

ggsave("../plots/joinplot.png", height = 10)
