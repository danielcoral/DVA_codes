## Own version of gmirror

library(tidyverse)
library(ggh4x)

mix <- rio::import("~/dva/files/mix.tsv")

bmi_t2d <- rio::import("~/dva/files/bmi_t2d.txt") %>%
    mutate(across(c(p.bmi, p.t2d),
                  function(x) -log10(ifelse(x < 1e-100, 1e-100, x))),
           chr = factor(chr, levels = 1:22)) %>%
    filter(p.bmi > -log10(0.05) | p.t2d > -log10(0.05)) %>%
    pivot_longer(c(p.bmi, p.t2d), names_to = "trait", values_to = "pval")

chroms <- bmi_t2d %>%
    group_by(chr) %>%
    summarise(l0 = min(hg19_pos) - 1,
              l1 = max(hg19_pos) + 1)

chrom_lengths <- chroms %>%
    transmute(chr, l = (l1 - l0) / 1e6)

chrom_lims <- chroms %>%
    pivot_longer(-chr, names_to = NULL, values_to = "hg19_pos")

## First trying with one
g <- bmi_t2d  %>%
    mutate(trait = recode(trait, "p.bmi" = "BMI", "p.t2d" = "T2D"),
           sig = ifelse(rsid %in% mix$rsid, hg19_pos, NaN)) %>%
    bind_rows(mutate(chrom_lims, trait = "BMI"),
              mutate(chrom_lims, trait = "T2D")) %>%
    group_by(trait) %>%
    group_map(
        ~{
            ggplot(.x, aes(pval, hg19_pos)) +
                geom_hline(aes(yintercept = sig), color = "grey",
                           alpha = .2) +
                geom_point(size = .1) +
                scale_x_continuous(breaks = c(25, 50, 75),
                                   expand = c(0,0)) +
                facet_wrap(~ chr, scales = "free_y", ncol = 1,
                           strip.position = "left") +
                labs(title = .y,
                     y = "Chromosome") +
                geom_vline(xintercept = -log10(5e-8),
                           linetype = "dashed", color = "red") +
                theme(axis.text.y = element_blank(),
                      axis.ticks.y = element_blank(),
                      axis.title.x = element_blank(),
                      axis.title.y = element_text(angle = 90),
                      axis.text.x = element_text(size = 6),
                      panel.spacing = unit(.1, "lines"),
                      panel.background = element_rect(fill = NA, color = "black"),
                      plot.title = element_text(hjust = .5),
                      strip.text.y.left = element_text(angle = 0, color = "black",
                                                       size = 6),
                      strip.background = element_blank(),
                      plot.margin = unit(c(0,0,0,0), "cm")) +
                force_panelsizes(rows = chrom_lengths$l)
        }
    ) %>%
    modify_at(2, ~ .x + theme(strip.text.y.left = element_blank(),
                              strip.background = element_blank(),
                              axis.title.y = element_blank())) %>%
    patchwork::wrap_plots(nrow = 1)

s <- mix %>%
    transmute(chr = factor(chr, levels = 1:22),
              hg19_pos, x = .5) %>%
    bind_rows(chrom_lims) %>%
    ggplot(aes(x, hg19_pos)) +
    geom_point(shape = 23, fill = "green", color = "black", size = .5) +
    facet_wrap(~chr, scales = "free_y", ncol = 1) +
    xlim(0,1) +
    labs(title = "Hits") +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank(),
          panel.spacing = unit(.1, "lines"),
          panel.background = element_rect(fill = NA, color = NA)) +
    force_panelsizes(rows = chrom_lengths$l)

gmirror <- list(g, s, grid::textGrob("-log10 p"), patchwork::plot_spacer()) %>%
    patchwork::wrap_plots(ncol = 2, nrow = 2, widths = c(.95, .05), heights = c(.99,.01))

save(gmirror, file = "../plot_files/gmirror.RData") 

ggsave("../docs/plots/gmirror.png", gmirror)
