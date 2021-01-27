## PRS and MR plot

library(tidyverse)
library(patchwork)

## Effects found
ukb_prs <- rio::import("../files/ukb_prs_res.tsv")
mr_res <- rio::import("../files/mr_res.tsv")

mrdat <- mr_res %>%
    group_by(disc, exposure, method) %>%
    mutate(sig = na_if(method == "ivw" & pval < 0.05/8, F) * 1,
           sigb0 = na_if(method == "intercept" & pval < 0.05, F) * 1,
           anym = na_if(any(duplicated(psi)), F) * 1) %>%
    filter(is.na(psi) | psi == 1.5) %>%
    group_by(method) %>%
    mutate(across(c(sig, sigb0, anym), ~.x * max(ci_up) * 1.2))
    

prs_mr_dat <- ukb_prs %>%
    filter(!trait %in% c("bmi", "t2d", "obesity", "hbp", "chd", "t_chol")) %>%
    transmute(trait, profile = gsub("prs_", "", term),
              b = estimate, p = p.value, ci_lo = lo, ci_up = up, method = "prs") %>%
    bind_rows(mrdat %>%
              transmute(trait = gsub("mass", "", exposure),
                        profile = ifelse(disc == 1, "disc", "conc"),
                        b, ci_lo, ci_up, p = pval, method,
                        sig, sigb0, anym)) %>%
    filter(method != "egger") %>%
    mutate(trait = factor(
               str_replace_all(trait,
                               c("sbp" = "SBP", "dbp" = "DBP",
                                 "hdl" = "HDL", "ldl" = "LDL", "tg" = "TG",
                                 "armfat" = "Arm fat mass", "legfat" = "Leg fat mass",
                                 "whr" = "WHR")),
               levels = c("SBP", "DBP", "HDL", "LDL", "TG",
                          "Arm fat mass", "Leg fat mass", "WHR")),
           method = recode_factor(method, prs = "PRS",
                                  ivw = "IVW \u03b2",
                                  conmix = "ConMix \u03b2",
                                  intercept = "Egger \u03b2\u2080"))

plotlist <- prs_mr_dat %>%
    group_by(method) %>%
    group_map(~{
        p <- .x %>%
            ggplot(aes(b, profile)) +
            geom_errorbar(aes(xmin = ci_lo, xmax = ci_up),
                          alpha = .4, size = .4, width = .2) +
            geom_point(aes(colour = profile), alpha = .5) +
            geom_point(shape = 1, col = "black", alpha = .4) +
            scale_colour_manual(values = c("red", "blue"),
                                labels = c("Concordant", "Discordant"),
                                guide = NULL) +
            geom_vline(xintercept = 0, linetype = "dashed") +
            facet_wrap(trait~., strip.position = "left", ncol = 1) +
            labs(title = .y$method) +
            theme(plot.title = element_text(hjust = .5),
                  axis.title.y = element_blank(),
                  axis.title.x = element_text(size = 8),
                  axis.text.x = element_text(size = 6),
                  axis.text.y = element_blank(),
                  axis.ticks.y = element_blank())
        if (.y == "PRS") {
            p + theme(strip.text.y.left = element_text(angle = 0))
        } else {
            p + theme(strip.background = element_blank(),
                      strip.text.y.left = element_blank()) +
                geom_point(aes(sig, profile), shape = 8, size = .8) +
                geom_point(aes(sigb0, profile), shape = "\u2605", size = 2) +
                geom_point(aes(anym, profile), shape = 18, size = 1.8)
        }
    })

p1 <- plotlist[[1]] + xlab("SD per 1kg/m\u00B2")
p2 <- wrap_plots(plotlist[2:3]) & xlab("logOR per SD")
p3 <- plotlist[[4]] + xlab("") + scale_x_continuous(breaks = c(-0.05, 0, 0.05))

prs_mr_plot <- wrap_plots(
    p1, p2, p3,
    widths = c(1,2,.5)
)
