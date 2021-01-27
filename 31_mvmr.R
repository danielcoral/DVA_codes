## Multivariable MR

library(tidyverse)

## Insulin DI and sensitivity index (Prokopenko et al)
ins_files <- paste0("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/",
                    c(paste0("ProkopenkoI_24699409_GCST005222/harmonised/",
                             "24699409-GCST005222-EFO_0006832.h.tsv.gz"),
                      paste0("ProkopenkoI_24699409_GCST005221/harmonised/",
                             "24699409-GCST005221-EFO_0004471.h.tsv.gz")))

## Adiponectin (Dastani et al)
adip_url <- paste0("https://www.mcgill.ca/genepi/files/genepi/",
                   "adipogen.discovery.eur_.meta_.public.release.part")
adip_files <- c(paste0(adip_url, 1:5, "_.txt"),
                paste0(adip_url, 6, "__0.txt"))

## Download
destfolder <- "../files/mr_ss/"

download.file(ins_files, method = "libcurl",
              destfile = paste0(destfolder, c("di", "isi"), ".tsv.gz"))

if(!dir.exists(paste0(destfolder, "adip")))
    dir.create(paste0(destfolder, "adip"))

download.file(adip_files, method = "libcurl",
              destfile = paste0(destfolder, "adip/", "adip", 1:6, ".txt"))

## Importing data
di <- vroom::vroom("../files/mr_ss/di.tsv.gz") %>%
    select(rsid = variant_id, ea_di = effect_allele,
           nea_di = other_allele, eaf_di = effect_allele_frequency,
           b_di = beta, se_di = standard_error, p_di = p_value) %>%
    filter(nchar(ea_di) == 1, nchar(nea_di) == 1,
           eaf_di > 0.01, eaf_di < 0.99)

isi <- vroom::vroom("../files/mr_ss/isi.tsv.gz") %>%
    select(rsid = variant_id, ea_isi = effect_allele,
           nea_isi = other_allele, eaf_isi = effect_allele_frequency,
           b_isi = beta, se_isi = standard_error, p_isi = p_value) %>%
    filter(nchar(ea_isi) == 1, nchar(nea_isi) == 1,
           eaf_isi > 0.01, eaf_isi < 0.99)

adip <- vroom::vroom(paste0("../files/mr_ss/adip/adip", 1:6, ".txt")) %>%
    select(rsid = marker, maf_adip = maf,
           ea_adip = reference_allele, nea_adip = other_allele,
           b_adip = beta, se_adip = se, p_adip = pvalue) %>%
    filter(nchar(ea_adip) == 1, nchar(nea_adip) == 1, maf_adip > 0.01)
           

## Concordant and discordant instruments
sigbmi_t2d <- rio::import("../files/sigbmi_t2d.RData")
mr_ivs <- rio::import("../files/mr_ivs.tsv")

## First, testing DI, ISI and adiponectin separately
inner_join(di, sigbmi_t2d, by = "rsid") %>%
    filter(p_di < 5e-8) ## No valid instruments!

inner_join(isi, sigbmi_t2d, by = "rsid") %>%
    filter(p_isi < 5e-8 ) ## No valid instruments!

adip_mrdat <- inner_join(adip, sigbmi_t2d, by = "rsid") %>%
    filter(p_adip < 5e-8) ## Valid instruments found!

adip_mrdat_h <- adip_mrdat %>%
    mutate(harmon = case_when(
               ea == ea_adip & nea == nea_adip ~ 1,
               ea == nea_adip & nea == ea_adip ~ -1,
               T ~ 0),
           b_adip = harmon * b_adip) %>%
    filter(harmon != 0) %>%
    group_by(orig_rsid) %>%
    arrange(is_proxy, r2, p_adip) %>%
    slice(1) %>%
    ungroup()

## Here we only found valid instruments for adiponectin
## And only for the discordant set
## Performing MR for adiponectin
library(MendelianRandomization)

adip_mr <- adip_mrdat_h %>%
    group_by(disc) %>%
    group_map(~{
        to_clump <- transmute(.x, rsid, pval = p_adip)
        clump_dat <- ieugwasr::ld_clump(to_clump, clump_kb = 500, clump_r2 = 0.01,
                                        plink_bin = genetics.binaRies::get_plink_binary(),
                                        bfile = "../files/1kg_ref/EUR")
        dat <- filter(.x, rsid %in% clump_dat$rsid)
        i <- mr_input(bx = dat$b_adip, bxse = dat$se_adip,
                      by = dat$beta.t2d, byse = dat$se.t2d,
                      snps = dat$rsid, exposure = "adiponectin",
                      outcome = "T2D")
        list(i, mr_ivw(i), mr_egger(i),
             lapply(seq(1,2,.1),
                    function(psi)
                        mr_conmix(i, psi = psi *
                                         sd(i@betaY/i@betaX))))
    })

## The caveat is that these SNPs are not in cis with ADIPOQ
## Therefore, this result is not valid
## Saving for supplementary data
adip_mrdat_h %>%
    filter(rsid %in% adip_input[[1]][[1]]@snps) %>%
    select(chr, pos, rsid, ea, nea,
           beta.bmi, se.bmi, p.bmi,
           b_adip, se_adip, p_adip,
           beta.t2d, se.t2d, p.t2d) %>%
    rio::export("../files/adiponectin_instrum.xlsx")

## MVMR
di <- rename_with(di, ~gsub("_di", "2", .x)) ## %>%
## inner_join(rename_with(adip, ~gsub("_adip", "3", .x)), by = "rsid")

isi <- rename_with(isi, ~gsub("_isi", "2", .x))

mvmr_res <- data.frame(trait1 = c("sbp", "dbp", "whr"),
                       trait2 = c("di", "di", "isi"),
                       profile = c(1, 1, 0)) %>%
    rowwise %>%
    group_map(~{
        dat1 <- mr_ivs %>%
            filter(trait_exp == .x$trait1, pval.exposure < 5e-8, disc == .x$profile) %>%
            transmute(trait1 = trait_exp, orig_rsid, rsid, is_proxy, r2, ea, nea,
                      b1 = beta.exposure, se1 = se.exposure, p1 = pval.exposure,
                      b_out = beta.t2d, se_out = se.t2d)
        dat2 <- get(.x$trait2)
        traitnames <- unlist(str_split(.x$trait2, "_"))
        joindat <- inner_join(dat1, dat2, by = "rsid") %>%
            mutate(harmon2 = case_when(ea == ea2 & nea == nea2 ~ 1,
                                       ea == nea2 & nea == ea2 ~ -1,
                                       T ~ 0), b2 = harmon2 * b2)
        if(length(traitnames) == 2){
            joindat <- joindat %>%
                mutate(harmon3 = case_when(ea == ea3 & nea == nea3 ~ 1,
                                           ea == nea3 & nea == ea3 ~ -1,
                                           T ~ 0), b3 = harmon3 * b3,
                       harmon2 = harmon2 * harmon3)
        }
        joindat <- joindat %>%
            filter(harmon2 != 0) %>%
            group_by(orig_rsid) %>%
            arrange(is_proxy, r2, p1) %>%
            slice(1) %>%
            ungroup
        to_clump <- transmute(joindat, rsid, pval = p1)
        clump_dat <- ieugwasr::ld_clump(to_clump, clump_kb = 500, clump_r2 = 0.01,
                                        plink_bin = genetics.binaRies::get_plink_binary(),
                                        bfile = "../files/mrtot2d_ref/ref")
        mrdat <- filter(joindat, rsid %in% clump_dat$rsid)
        bx = as.matrix(select(mrdat, starts_with("b") & !ends_with("_out")))
        bxse = as.matrix(select(mrdat, starts_with("se") & !ends_with("_out")))
        i <- mr_mvinput(bx = bx, bxse = bxse, by = mrdat$b_out, byse = mrdat$se_out)
        resivw <- mr_mvivw(i); resegger <- mr_mvegger(i); resmed <- mr_mvmedian(i)
        bind_rows(tibble(method = "ivw", trait = c(.x$trait1, traitnames),
                         b = resivw@Estimate, se = resivw@StdError,
                         ci_lo = resivw@CILower, ci_up = resivw@CIUpper,
                         p = resivw@Pvalue, Q = resivw@Heter.Stat[1],
                         Qdf = resivw@SNPs - length(resivw@Exposure),
                         Qp = resivw@Heter.Stat[2]),
                  tibble(method = "median", trait = c(.x$trait1, traitnames),
                         b = resmed@Estimate, se = resmed@StdError,
                         ci_lo = resmed@CILower, ci_up = resmed@CIUpper,
                         p = resmed@Pvalue),
                  tibble(method = "egger", trait = c("intercept", .x$trait1, traitnames),
                         b = c(resegger@Intercept, resegger@Estimate),
                         se = c(resegger@StdError.Int, resegger@StdError.Est),
                         ci_lo = c(resegger@CILower.Int, resegger@CILower.Est),
                         ci_up = c(resegger@CIUpper.Int, resegger@CIUpper.Est),
                         p = c(resegger@Pvalue.Int, resegger@Pvalue.Est),
                         Q = resegger@Heter.Stat[1],
                         Qdf = resegger@SNPs - length(resegger@Exposure) - 1,
                         Qp = resegger@Heter.Stat[2])) %>%
            mutate(test = paste(.x$trait1, "adj"), profile = .x$profile)
    }) %>%
    bind_rows

rio::export(mvmr_res, "../files/mvmr_res.tsv")

mvmr_res %>%
    rename_with(toupper) %>%
    rename(CI_LOWER = CI_LO, CI_UPPER = CI_UP,
           BETA = B, MODEL = TEST) %>%
    mutate(METHOD = gsub("Ivw", "IVW", str_to_sentence(METHOD)),
           TRAIT = str_replace_all(TRAIT,c("sbp" = "SBP", "dbp" = "DBP",
                                           "di" = "Insulin disposition index",
                                           "intercept" = "Intercept", "whr" = "WHR",
                                           "isi" = "Insulin sensitivity index")),
           MODEL = paste0(toupper(MODEL), "USTED"),
           PROFILE = ifelse(PROFILE == 1, "Discordant", "Concordant"),
           across(c(BETA, SE, CI_LOWER, CI_UPPER), exp)) %>%
    relocate(MODEL, PROFILE) %>%
    rio::export("../files/mvmr_res.xlsx")

mvmr_res <- rio::import("../files/mvmr_res.tsv")

mvmr_plot <- list(
    mvmr_res %>%
    filter(trait != "intercept") %>%
    mutate(method = recode_factor(method,
                                  egger = "Egger", median = "Median", ivw = "IVW"),
           trait = recode_factor(trait,
                                 sbp = "SBP", dbp = "DBP", di = "Disposition Index", whr = "WHR",
                                 isi = "Sensitivity Index"),
           test = recode_factor(test,
                                `sbp adj` = "SBP MR",
                                `dbp adj` = "DBP MR",
                                `whr adj` = "WHR MR")) %>%
    ggplot(aes(b, method)) +
    geom_point() +
    geom_errorbar(aes(xmin = ci_lo, xmax = ci_up),
                  alpha = .4, size = .4, width = .2) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = .2) +
    facet_wrap(test ~ trait, dir = "v", ncol = 3) +
    labs(x = "logOR per SD", y = NULL),
    mvmr_res %>%
    filter(trait == "intercept") %>%
    mutate(test = recode_factor(test,
                                `sbp adj` = "SBP MR\nEgger Intercept",
                                `dbp adj` = "DBP MR\nEgger Intercept",
                                `whr adj` = "WHR MR\nEgger Intercept")) %>%
    ggplot(aes(b, method)) +
    geom_point() +
    geom_errorbar(aes(xmin = ci_lo, xmax = ci_up),
                  alpha = .4, size = .4, width = .2) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = .2) +
    facet_wrap(~test) +
    theme(axis.title = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))
) %>%
    patchwork::wrap_plots(ncol = 1, heights = c(9,1))

ggsave("../plots/mvmr_plot.png", mvmr_plot)
