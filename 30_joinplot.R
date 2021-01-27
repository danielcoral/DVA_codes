## Joining MR and PRS plots

unlink(list.files(pattern = ".~"))
lapply(list.files(pattern = "_plot_"), source)

library(cowplot)

joinplot <- plot_grid(plot_grid(prs_mr_plot, survival_plot, ncol = 2,
                                rel_widths = c(.8, .2), labels = c("A", "B")),
                      plot_grid(pw_biovu_plot, lw_biovu_plot,
                                ncol = 2, rel_widths = c(.6,.4), labels = c("C", "D")),
                      ncol = 1, rel_heights = c(.3,.7))

ggsave("../plots/joinplot.png", height = 10)
