## BioVU PheWAS results

library(tidyverse)
library(ggh4x)

pw_biovu <- vroom::vroom(list.files("~/dva/files/biovu/phewas", pattern = "w_unit", full.names = T),
                   id = "profile") %>%
    drop_na %>%
    ## Threshold from Verma et al. 2018
    filter(n_cases > 200) %>%
    mutate(profile = gsub(".*phewas_(.*)ordant_GRS.*", "\\1_prs", profile),
           padj = p.adjust(p, "fdr")) %>%
    { ## Significant signals (FDR < 5%)
        sig_dx <- filter(., padj < 0.05) %>%
            pull(phecode) %>%
            unique
        filter(., phecode %in% sig_dx) } %>%
    ## Retaining synonyms of obesity
    mutate(obesityrel = description %in% c("Overweight, obesity and other hyperalimentation",
                                           "Obesity", "Morbid obesity", "Bariatric surgery")) %>%
    ## Calculating effect per BMI unit
    ##left_join(prs_in_bmi, by = c("profile" = "term")) %>%
    mutate(##prs_sd = log(OR.per.sd) / beta,
           ##beta = beta * prs_sd / u,
           ##SE = SE * prs_sd / u,
           conf.high = beta + qnorm(1-0.05/2)*SE,
           conf.low = beta - qnorm(1-0.05/2)*SE) %>%
    ## Comparing effects of profiles
    group_by(phecode) %>%
    mutate(bcomp = abs(beta[1] - beta[2]),
           secomp = sqrt(SE[1]^2 + SE[2]^2),
           pcomp = pnorm(-bcomp/secomp)) %>%
    ungroup %>%
    ## Traits with significant difference (FDR < 5%)
    mutate(padj = p.adjust(pcomp, "fdr")) %>%
    filter(obesityrel | padj < 0.05) %>%
    ## What happens with these traits after excluding T1D patients?
    mutate(dataset = "all") %>%
    {
        tokeep <- pull(., phecode)
        not1d <- vroom::vroom(list.files("~/dva/files/biovu/phewas_rm_t1d",
                                pattern = "w_unit", full.names = T),
                     id = "profile") %>%
            filter(phecode %in% tokeep) %>%
            mutate(profile = gsub(".*phewas_(.*)ordant_GRS.*", "\\1_prs", profile)) %>%
            ## Calculating effect per BMI unit
            ##left_join(prs_in_bmi, by = c("profile" = "term")) %>%
            mutate(##prs_sd = log(OR.per.sd) / beta,
                   ##beta = beta * prs_sd / u,
                   ##SE = SE * prs_sd / u,
                   conf.high = beta + qnorm(1-0.05/2)*SE,
                   conf.low = beta - qnorm(1-0.05/2)*SE) %>%
            ## Comparing effects of profiles
            group_by(phecode) %>%
            mutate(bcomp = abs(beta[1] - beta[2]),
                   secomp = sqrt(SE[1]^2 + SE[2]^2),
                   pcomp = pnorm(-bcomp/secomp)) %>%
            ungroup %>%
            ## Traits that continue with significant difference (FDR < 5%)
            mutate(dataset = "not1d", padj = p.adjust(pcomp, "fdr")) %>%
            mutate(obesityrel = description %in% c("Overweight, obesity and other hyperalimentation",
                                                   "Obesity", "Morbid obesity", "Bariatric surgery"))
        bind_rows(., not1d) }

rio::export(pw_biovu, "~/dva/files/pw_biovu.tsv")

## Editing before plotting
pw_biovu_dat <- pw_biovu %>%
    mutate(dataset = ifelse(dataset == "all", "All", "No T1D"),
           ## Signals that survived exclusion of T1D patients
           sigdif = na_if(dataset == "No T1D" & padj < 0.05, FALSE),
           sigdif = sigdif * max(conf.high) * 1.05,
           ## Shortening names
           description = str_replace_all(
               description,
               c("Hypertension|Hypertensive" = "HT",
                 "chronic kidney disease|Chronic Kidney Disease" = "CKD",
                 "Chronic renal failure \\[CKD\\]" = "CKD",
                 "Peripheral vascular disease|peripheral circulatory disorders" = "PAD",
                 "diseases|disease|disorders|disorder|manifestations" = "dis",
                 "Secondary diabetes mellitus" = "Secondary diabetes",
                 "Sensorineural hearing loss" = "Hearing loss",
                 "; interstitial and compensatory emphysema" = "",
                 " and/or " = "/",
                 " \\(intermediate coronary syndrome\\)" = "",
                 ", unspecified" = "",
                 ", except cervix" = "",
                 "or abnormal glucose tolerance complicating" = "in",
                 " of the extremities" = " - limbs",
                 "Overweight.*" = "Overweight",
                 "Type 2 diabetes|Diabetes type 2" = "T2D",
                 "\\(CHF\\) " = "",
                 "Congestive heart failure" = "CHF",
                 "Type 1 diabetes" = "T1D",
                 " with " = " - ",
                 "Polyneuropathy in diabetes" = "Diabetic neuropathy",
                 "Disorders of lipoid metabolism" = "Lipoid metabolism dis",
                 "Hyperlipidemia" = "Hyperlipidaemia",
                 "Ischemic" = "Ischaemic",
                 "Staphylococcus aureus" = "Staph aureus",
                 "Complications of cardiac.*" = "Cardiac/vascular device complications")
           ),
           ## Grouping
           priority = factor(case_when(obesityrel ~ "Obesity",
                                       grepl("t1d|t2d|diabet|glucose|insulin",
                                             description, ignore.case = TRUE) ~ "Diabetes",
                                       group == "genitourinary" | grepl("CKD", description) ~ "Renal",
                                       group == "circulatory system" ~ "CV",
                                       grepl("lipid|lipoid|cholesterol",
                                             description, ignore.case = TRUE) ~ "Lipids",
                                       TRUE ~ "Other"),
                             levels = c("Obesity", "Diabetes", "CV",
                                        "Lipids", "Renal", "Other")))

## Plot
pwbiovu_plot <- pw_biovu_dat %>%
    ggplot(aes(beta,
               interaction(description,
                           reorder(priority, desc(priority))),
               group = profile)) +
    geom_errorbar(aes(xmin = conf.low, xmax = conf.high),
                  alpha = .4, size = .4, width = .2,
                  position = position_dodge(width = 0.5)) +
    geom_point(aes(color = profile), alpha = .5,
               position = position_dodge(width = 0.5)) +
    geom_point(shape = 1, col = "black", alpha = .4,
               position = position_dodge(width = 0.5)) +
    geom_point(aes(x = sigdif), shape = "*", size = 4) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    ##    scale_color_manual(values = c("red", "blue"), guide = "none") +
    scale_colour_discrete(labels = c("Concordant", "Discordant"),
                          guide = guide_legend(title = NULL)) +
    guides(y = "axis_nested") +
    facet_wrap(~ dataset) +
    labs(x = "logOR per allele") +
    theme(axis.title.y = element_blank(),
          legend.position = "top",
          axis.text.y = element_text(size = 7),
          ggh4x.axis.nesttext.y = element_text(angle = 90, hjust = .5,
                                               margin = margin(l = 10),
                                               face = "bold"),
          strip.text = element_text(face = "bold"))

ggsave("../plots/biovu_pw.png")

## Saving to include in supplementary data
pw_biovu_dat %>%
    transmute(INCLUDED_T1D = recode(dataset, All = "YES",
                                    `No T1D` = "NO"),
              PRS = recode(profile, prs_conc = "Concordant",
                           prs_disc = "Discordant"),
              PHECODE = as.character(phecode), DESCRIPTION = description, P = p,
              N_CASES = n_cases, N_CONTROLS = n_controls,
              ## LOG OR
              BETA = beta, CI_LOWER = conf.low, CI_UPPER = conf.high,
              ## OR
              OR = exp(beta), CI_LOWER_OR = exp(conf.low), CI_UPPER_OR = exp(conf.high),
              P_DIFF = pcomp) %>%
    arrange(INCLUDED_T1D, PHECODE) %>%
    rio::export("~/dva/files/pw_biovu.xlsx")
