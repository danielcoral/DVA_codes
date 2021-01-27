## LabWAS plot

library(tidyverse)

lw_biovu <- rio::import("../files/lw_biovu.tsv")

lw_biovu_plot <- lw_biovu %>%
    mutate(trait = str_replace_all(
               Full_name,
               c("Hemoglobin A1c.*" = "HbA1c",
                 "Glucose lab" = "Glucose mmol/L",
                 "Glucose.*strip" = "Glucose - Strip",
                 "Glucose \\[.*?\\] in Blood" = "Glucose mg/dL",
                 ".*HDL.*"= "HDL",
                 ".*LDL.*" = "LDL",
                 "Cholesterol.*" = "Cholesterol",
                 "Urea nitrogen.*" = "BUN",
                 "^Erythrocyte distribution width (\\[.*?\\]) by.*" = "EDW \\1",
                 "Triglyceride.*" = "TG",
                 "MCV.*" = "MCV",
                 "^Creatinine.*" = "Creatinine",
                 "MCHC .*" = "MCHC",
                 "MCH .*" = "MCH",
                 "Calcium.*" = "Calcium",
                 "Anion gap.*" = "Anion gap",
                 "C reactive protein.*" = "CRP",
                 "Parathyrin.*" = "PTH",
                 "Platelet.*" = "Platelet volume",
                 "Leukocytes.*" = "Leucocyte count",
                 "Entitic " = "")
           ),
           priority = ifelse(pop == "AA", 0, p),
           trait = ifelse(pop == "AA", paste0(trait, "-AA"), trait)) %>%
    group_by(pop, trait) %>%
    mutate(bcomp = abs(beta[1] - beta[2]),
           secomp = sqrt(SE[1]^2 + SE[2]^2),
           pcomp = pnorm(-abs(bcomp/secomp))) %>%
    { ## Saving for supplementary data
        df <- ungroup(.)
        left_join(df, rio::import("../files/prs_in_bmi.tsv"),
                  by = c("profile" = "term")) %>%
            transmute(POPULATION = pop,
                      PRS = recode(profile, prs_conc = "Concordant",
                                   prs_disc = "Discordant"),
                      FULL_NAME = Full_name, P = p,
                      BETA = beta * u, CI_LOWER = ci_lo * u, CI_UPPER = ci_up * u,
                      BETA2 = beta, CI_LOWER2 = ci_lo, CI_UPPER2 = ci_up,
                      P_DIFF = pcomp) %>%
            arrange(POPULATION, FULL_NAME) %>%
            rio::export("../files/lw_biovu.xlsx")
        df } %>%
    group_by(pop) %>%
    filter(pcomp < 0.05/length(unique(trait))) %>%
    ggplot(aes(beta, profile)) +
    geom_errorbar(aes(xmin = ci_lo, xmax = ci_up),
                  alpha = .4, size = .4, width = .2) +
    geom_point(aes(colour = profile), alpha = .5) +    
    geom_point(shape = 1, col = "black", alpha = .4) +
    scale_colour_manual(values = c("red", "blue"),
                        labels = c("Concordant", "Discordant"),
                        guide = "none") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    facet_grid(reorder(trait, priority) ~ ., switch = "y") +
    labs(title = "", x = "SD per 1kg/m\u00B2") +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size = 8),
          axis.text.x = element_text(size = 6),
          axis.ticks.y = element_blank(),
          strip.text.y.left = element_text(angle = 0))
