## Plotting effects of concordant and discordant PRS in BioVU

library(tidyverse)

pw_biovu <- rio::import("../files/pw_biovu.tsv")

pw_biovu_plot <- pw_biovu %>%
    mutate(description = str_replace_all(
               description,
               c("Overweight.*" = "Overweight",
                 "Type 2 diabetes" = "T2D",
                 "\\(CHF\\) " = "",
                 "Congestive heart failure" = "CHF",
                 "Type 1 diabetes" = "T1D",
                 " with " = " - ",
                 "manifestations" = "",
                 "Polyneuropathy in diabetes" = "Diabetic neuropathy",
                 "Disorders of lipoid metabolism" = "Disorders lipoid metab",
                 "Hyperlipidemia" = "Hyperlipidaemia",
                 "Ischemic" = "Ischaemic")
           ),
           priority = case_when(
               description %in% c("Obesity", "T2D") ~ 1,
               description %in% c("Morbid obesity", "Overweight",
                                  "Bariatric surgery") ~ 2,
               grepl("diabet|t1d|t2d|insulin", description,
                     ignore.case = T) ~ 3,
               grepl("hypertension", description,
                     ignore.case = T) ~ 4,
               grepl("lipid|lipoid", description,
                     ignore.case = T) ~ 5,
               grepl("coronary|heart|chf|cardiomegaly|myocardial",
                     description, ignore.case = T) ~ 6,
               grepl("renal", description, ignore.case = T) ~ 7,
               grepl("bacter|infect", description, ignore.case = T) ~ 8,
               T ~ 9
           ),
           dataset = ifelse(dataset == "all", "All", "No T1D")) %>%
    group_by(dataset, phecode) %>%
    mutate(bcomp = abs(beta[1] - beta[2]),
           secomp = sqrt(SE[1]^2 + SE[2]^2),
           pcomp = pnorm(-abs(bcomp/secomp))) %>%
    group_by(dataset) %>%
    mutate(sigdif = na_if(pcomp < 0.05/length(unique(phecode)), F) *
               max(ci_up) * 1.05) %>%
    group_by(dataset, phecode) %>%
    mutate(sigdif = ifelse(duplicated(phecode), NA, sigdif)) %>%
    { ## Saving to include in supplementary data
        df <- .
        left_join(ungroup(df), rio::import("../files/prs_in_bmi.tsv"),
                  by = c("profile" = "term")) %>%
            transmute(INCLUDED_T1D = recode(dataset, All = "YES",
                                            `No T1D` = "NO"),
                      PRS = recode(profile, prs_conc = "Concordant",
                                   prs_disc = "Discordant"),
                      PHECODE = as.character(phecode), DESCRIPTION = description, P = p,
                      ## LOG OR PER SD
                      BETA = beta * u, CI_LOWER = ci_lo * u, CI_UPPER = ci_up * u,
                      ## LOG OR PER BMI UNIT
                      BETA2 = beta, CI_LOWER2 = ci_lo, CI_UPPER2 = ci_up,
                      ## OR PER BMI UNIT
                      BETA3 = exp(beta), CI_LOWER3 = exp(ci_lo), CI_UPPER3 = exp(ci_up),
                      P_DIFF = pcomp) %>%
            arrange(INCLUDED_T1D, PHECODE) %>%
            rio::export("../files/pw_biovu.xlsx")
        df } %>%
    ggplot(aes(beta, profile)) +
    geom_errorbar(aes(xmin = ci_lo, xmax = ci_up),
                  alpha = .4, size = .4, width = .2) +
    geom_point(aes(color = profile), alpha = .5) +
    geom_point(shape = 1, col = "black", alpha = .4) +
    geom_point(aes(sigdif, profile), shape = 8, size = .8,
               position = position_nudge(y = 0.5)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_color_manual(values = c("red", "blue"),
                       guide = "none") +
    facet_grid(reorder(description, priority) ~ dataset, switch = "y") +
    labs(x = "logOR per 1kg/m\u00B2") +
    theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          strip.text.y.left = element_text(angle = 0),
          axis.title.x = element_text(size = 8),
          axis.text.x = element_text(size = 6))
