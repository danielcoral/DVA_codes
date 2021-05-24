## LabWAS plot

library(tidyverse)

## LabWAS
lw_biovu <- vroom::vroom(list.files("~/dva/files/biovu/labwas", full.names = T),
                         id = "path") %>%
    filter(Predictor == "Unit_PRS") %>% ## Unweighted PRS
    ## Significant signals (FDR < 5%)
    mutate(pop_profile = gsub(".*labwas/(.*)_PRS.*", "\\1", path),
           path = NULL,
           padj = p.adjust(p, "fdr")) %>%
    separate(pop_profile, c("pop", "profile")) %>%
    {
        sig_lab_EUR <- filter(., pop == "EA" & padj < 0.05) %>%
            pull(Lab) %>%
            unique
        sig_lab_AFR <- filter(., pop == "AA" & padj < 0.05) %>%
            pull(Lab) %>%
            unique
        filter(.,
               pop == "EA" & Lab %in% sig_lab_EUR |
               pop == "AA" & Lab %in% sig_lab_AFR) } %>%
    ## Calculating effect per BMI unit
    mutate(profile = tolower(paste0(gsub("ordant", "", profile), "_prs"))) %>%
    ##left_join(prs_in_bmi, by = c("profile" = "term")) %>%
    mutate(##beta = Estimate / u,  ##SE = SE / u,
        beta = Estimate,
        conf.high = beta + qnorm(1-0.05/2)*SE,
        conf.low = beta - qnorm(1-0.05/2)*SE) %>%
    ## Comparing effects of profiles
    group_by(pop, Lab) %>%
    mutate(bcomp = abs(beta[1] - beta[2]),
           secomp = sqrt(SE[1]^2 + SE[2]^2),
           pcomp = pnorm(-bcomp/secomp)) %>%
    ungroup %>%
    ## Traits with significant difference (FDR < 5%)
    mutate(padj = p.adjust(pcomp, "fdr")) %>%
    filter(padj < 0.05)

rio::export(lw_biovu, "~/dva/files/lw_biovu.tsv")

lw_biovu_dat <- lw_biovu %>%
    mutate(trait = recode_factor(
               Full_name,
               `Hemoglobin A1c (Glycated) Hemoglobin A1c (Glycated)` = "HbA1c",
               `Glucose lab` = "Glucose (mmol/L)",
               `Glucose [Mass/volume] in Blood` = "Glucose (mg/dL)",
               `Glucose [Mass/volume] in Blood by Automated test strip` = "Glucose (trip)",
               `Glucose [Mass/volume] in Serum or Plasma --2 hours post 100 g glucose PO` = "2h OGTT",
               `Anion gap serum/plasma` = "Anion gap",
               `Cholesterol [Mass/volume] in Serum or Plasma` = "Cholesterol",
               `Cholesterol in HDL [Mass/volume] in Serum or Plasma` = "HDL",
               `Cholesterol in LDL [Mass/volume] in Serum or Plasma by calculation` = "LDL (calc)",
               `Triglyceride [Mass/volume] in Serum or Plasma` = "TG",
               `Erythrocyte distribution width [Ratio] by Automated count Erythrocyte distribution width [Ratio] by Automated count` = "EDW ratio",
               `Erythrocyte distribution width [Entitic volume] by Automated count` = "EDW volume",
               `MCV [Entitic volume] by Automated count` = "MCV",
               `MCH [Entitic mass] by Automated count` = "MCH",
               `MCHC [Mass/volume] by Automated count` = "MCHC",
               `Erythrocyte sedimentation rate by Westergren method` = "ESR",
               `Leukocytes [#/volume] in Blood by Automated count` = "#Leuco",
               `Platelet mean volume [Entitic volume] in Blood by Automated count Platelet mean volume [Entitic volume] in Blood by Automated count` = "Platelet volume",
               `Urea nitrogen [Mass/volume] in Urine` = "Urea (urine)",
               `Urea nitrogen serum/plasma` = "Urea (blood)",
               `Creatinine [Mass/volume] in Blood Creatinine [Mass/volume] in Blood` = "Creatinine",
               `Phosphate` = "Phosphate",
               `C reactive protein [Mass/volume] in Serum or Plasma` = "CRP",
               `Parathyrin.intact [Mass/volume] in Serum or Plasma` = "PTH"
           ))

lwbiovu_plot <- lw_biovu_dat %>%
    ggplot(aes(beta, reorder(trait, desc(trait)), group = profile)) +
    geom_errorbar(aes(xmin = conf.low, xmax = conf.high),
                  alpha = .4, size = .4, width = .2,
                  position = position_dodge(0.5)) +
    geom_point(aes(colour = profile), alpha = .5,
               position = position_dodge(0.5)) +    
    geom_point(shape = 1, col = "black", alpha = .4,
               position = position_dodge(0.5)) +
    ##    scale_colour_manual(values = c("red", "blue"), guide = "none") +
    scale_colour_discrete(labels = c("Concordant", "Discordant"),
                          guide = guide_legend(title = NULL)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    facet_grid(pop ~ ., scales = "free_y", space = "free_y") +
    labs(title = "", x = "SD per allele") +
    theme(axis.title.y = element_blank(),
          legend.position = "top",
          strip.text.y = element_text(angle = 0))

ggsave("../plots/lw_biovu.png")

## Saving for supplementary data
lw_biovu_dat %>%
    transmute(POPULATION = pop,
              PRS = recode(profile, prs_conc = "Concordant",
                           prs_disc = "Discordant"),
              FULL_NAME = Full_name, N, P = p,
              BETA = beta, CI_LOWER = conf.low, CI_UPPER = conf.high,
              P_DIFF = pcomp) %>%
    arrange(POPULATION, FULL_NAME) %>%
    rio::export("~/dva/files/lw_biovu.xlsx")
