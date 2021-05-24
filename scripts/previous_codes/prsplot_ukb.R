## PRS and MR plot

library(tidyverse)
library(ggh4x)

## Effects found
ukb_prs <- rio::import("../files/ukb_prs_res.tsv") 

prs_dat <- ukb_prs %>%
    filter(!trait %in% c("bmi0", "t2d")) %>%
    mutate(trait = recode_factor(trait,
                                 sbp0 = "SBP", dbp0 = "DBP",
                                 hdl0 = "HDL", ldl0 = "LDL", tg0 = "TG",
                                 armfat0 = "Arm fat mass", legfat0 = "Leg fat mass",
                                 whr0 = "WHR"),
           term = recode(term, conc_prs = "Concordant", disc_prs = "Discordant"),
           trait_group = factor(case_when(trait %in% c("SBP", "DBP") ~ "BP",
                                          trait %in% c("HDL", "LDL", "TG")  ~ "Lipids",
                                          trait %in% c("Arm fat mass", "Leg fat mass") ~ "Peripheral\nfat mass",
                                          trait == "WHR" ~ "Central adiposity")),
           sigdif = na_if(padj < 0.05, FALSE),
           sigdif = sigdif * max(conf.high) * 1.05)

prs_plot <- prs_dat %>%
    ggplot(aes(estimate,
               weave_factors(reorder(trait, desc(trait)), trait_group),
               group = term)) +
    geom_errorbar(aes(xmin = conf.low, xmax = conf.high),
                  alpha = .4, size = .4, width = .2,
                  position = position_dodge(width = 0.25)) +
    geom_point(aes(colour = term), alpha = .5,
               position = position_dodge(width = 0.25)) +
    geom_point(shape = 1, col = "black", alpha = .4,
               position = position_dodge(width = 0.25)) +
    geom_point(aes(x = sigdif), shape = "*", size = 4) +
    scale_colour_manual(values = c("red", "blue"), guide = "none") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(x = "SD per 1 kg/m\u00b2") +
    guides(y = "axis_nested") +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_text(hjust = .5),
          ggh4x.axis.nesttext.y = element_text(face = "bold"),
          strip.text = element_text(face = "bold"))
