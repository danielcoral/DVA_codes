## Supplementary data

library(tidyverse)

## List of concordant and discordant SNPs
mix <- rio::import("../files/mix.txt")

mix %>%
    mutate(across(c(p.bmi, p.t2d),
                  ~ paste0(round(10^(log10(.x) %% 1), 1), "x10^", log10(.x) %/% 1)),
           profile = ifelse(disc == 0, "Concordant", "Discordant"),
           disc = NULL) %>%
    arrange(desc(profile)) %>%
    rio::export("../files/mix.xlsx")

## List of traits in the SNP-Trait Matrix
premat_phen <- rio::import("../files/premat_phen.txt")

premat_phen %>%
    select(trait, study, pmid, ancestry, year,
           n, n_cases, n_controls, unit, dataset) %>%
    group_by(trait, dataset) %>%
    mutate(max_n = max(n), max_ncases = max(n_cases),
           max_ncontrols = max(n_controls)) %>%
    ungroup() %>%
    select(-c(n, n_cases, n_controls)) %>%
    unique() %>%
    arrange(trait) %>%
    rio::export("../files/premat_phen.xlsx")

## Complete results of PRS in UKB
std_devs <- rio::import("../files/ukb_f.RData") %>%
    summarise(across(all_of(ukb_prs %>%
                            filter(grepl("^sd", units)) %>%
                            pull(trait)),
                     sd, na.rm = T)) %>%
    pivot_longer(everything(), names_to = "trait", values_to = "sd_value")

prs_in_bmi <- rio::import("../files/prs_in_bmi.tsv")

ukb_prs <- rio::import("../files/ukb_prs_res.tsv") %>%
    left_join(prs_in_bmi, by = "term") %>%
    left_join(std_devs, by = "trait") %>%
    filter(!trait %in% c("t_chol", "hbp", "obesity", "chd")) %>%
    rename_with(~toupper(recode(.x, lo = "ci_lower", up = "ci_upper")),
                c(estimate, lo, up)) %>%
    transmute(
        PRS = recode(term, prs_conc = "Concordant", prs_disc = "Discordant"),
        TRAIT = factor(recode(toupper(trait), OBESITY = "Obesity",
                              ARMFAT = "Arm fat mass", LEGFAT = "Leg fat mass"),
                       levels = c("BMI", "T2D", "SBP", "DBP", "HDL", "TG", "LDL",
                                  "Arm fat mass", "Leg fat mass", "WHR")),
        UNITS = case_when(TRAIT == "BMI" ~ "KG/M2",
                          TRAIT == "T2D" ~ "LOGOR",
                          TRAIT %in% c("SBP", "DBP") ~ "MMHG",
                          TRAIT %in% c("HDL", "LDL", "TG") ~ "MMOL/L",
                          grepl("fat mass", TRAIT) ~ "KG",
                          TRAIT == "WHR" ~ "RATIO"),
        SD = sd_value,
        P = (log(2) +
             pnorm(abs(statistic), lower.tail = F, log.p = T)) /
            log(10),
        P = paste(round(10^(P %% 1), 2), P %/% 1, sep = "E"),
        across(c(ESTIMATE, CI_LOWER, CI_UPPER),
               ~ ifelse(trait != "BMI", .x * u, .x),
               .names = "{.col}_sd.prssd"),
        across(c(ESTIMATE, CI_LOWER, CI_UPPER),
               ~ ifelse(trait != "BMI", .x * sd_value * u, u),
               .names = "{.col}_nat.prssd"),
        across(c(ESTIMATE, CI_LOWER, CI_UPPER),
               ~ ifelse(trait != "BMI", .x, NaN),
               .names = "{.col}_sd.predbmi"),
        across(c(ESTIMATE, CI_LOWER, CI_UPPER),
               ~ ifelse(trait != "BMI", .x * sd_value, NaN),
               .names = "{.col}_nat.predbmi")
    ) %>%
    arrange(PRS, TRAIT)

rio::export(ukb_prs, "../files/ukb_prs_res.xlsx")

## MR
mr_res <- rio::import("../files/mr_res.tsv")

mr_res %>%
    filter(is.na(psi) | psi == 1.5) %>%
    transmute(disc = ifelse(disc == 0, "", "2"),
              EXPOSURE = toupper(exposure),
              METHOD = recode(toupper(method), INTERCEPT = "EGGER_B0"),
              N_SNPS = nsnp,
              BETA = exp(b), CI_LOWER = exp(ci_lo), CI_UPPER = exp(ci_up),
              P = pval) %>%
    pivot_wider(names_from = disc,
                values_from = c(N_SNPS, BETA, CI_LOWER, CI_UPPER, P),
                names_sep = "") %>%
    relocate(ends_with("2"), .after = last_col()) %>%
    rio::export("../files/mr_res.xlsx")

## ConMix
mr_res %>%
    filter(method == "conmix") %>%
    transmute(SUBTYPE = ifelse(disc == 0, "Concordant", "Discordant"),
              EXPOSURE = toupper(exposure),
              BETA = exp(b), P = pval,
              PSI_FACTOR = psi, PSI = psitrue,
              CI_LOWER = exp(ci_lo), CI_UPPER = exp(ci_up)) %>%
    arrange(PSI, desc(BETA > CI_LOWER & BETA < CI_UPPER)) %>%
    mutate(DISJOINT_CI = duplicated(PSI)) %>%
    rio::export("../files/conmix_res.xlsx")
