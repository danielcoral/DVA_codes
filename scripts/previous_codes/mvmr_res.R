## Multivariable MR

library(tidyverse)
library(MendelianRandomization)

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

### Are they within the ADIPOQ gene?
adipoq_chr <- 3
adipoq_pos <- c(186560463,186576252)

adip_mrdat %>%
    filter(chr == adipoq_chr,
           pos > adipoq_pos[1] - 500000,
           pos < adipoq_pos[2] - 500000) ## No!
## This means they are unsuitable for MR

## MVMR - Data
di <- rename_with(di, ~gsub("_di", "2", .x))
isi <- rename_with(isi, ~gsub("_isi", "2", .x))

## Custom MVMR Egger function for preventing change in orientation
mvmr.eggerfx <- function(Bx, By, Bxse, Byse) {
    summary <- summary(lm(By ~ Bx, weights = Byse^(-2)))
    data.frame(method = c("intercept",
                          rep("egger", dim(Bx)[2])),
               b = summary$coef[,1],
               se = summary$coef[,2]/min(summary$sigma,1),
               nsnp = dim(Bx)[1]) %>%
        mutate(pvalue = 2*pnorm(-abs(b/se)),
               Q_df = dim(Bx)[1] - dim(Bx)[2] - 1,
               Q = Q_df * (summary$sigma^2),
               Q_pval = pchisq(Q, df = Q_df, lower.tail = FALSE),
               ci_lo = b - qnorm(1-0.05/2) * se,
               ci_up = b + qnorm(1-0.05/2)* se)
}

## Results
mvmr_res <- map2(
    list("sbp", "dbp", "whr"),
    list("di", "di", "isi"),
    function(trait1, trait2) {
        dat1 <- mr_ivs %>%
            filter(trait_exp == trait1, pval.exposure < 5e-8) %>%
            transmute(trait1 = trait_exp, orig_rsid, rsid, is_proxy, r2, ea, nea,
                      b1 = beta.exposure, se1 = se.exposure, p1 = pval.exposure,
                      b_out = beta.t2d, se_out = se.t2d)
        dat2 <- get(trait2)
        joindat <- inner_join(dat1, dat2, by = "rsid") %>%
            mutate(harmon2 = case_when(ea == ea2 & nea == nea2 ~ 1,
                                       ea == nea2 & nea == ea2 ~ -1,
                                       T ~ 0), b2 = harmon2 * b2) %>%            
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
        bx <- as.matrix(select(mrdat, starts_with("b") & !ends_with("_out")))
        bxse <- as.matrix(select(mrdat, starts_with("se") & !ends_with("_out")))
        i <- mr_mvinput(bx = bx, bxse = bxse, by = mrdat$b_out, byse = mrdat$se_out)
        resivw <- mr_mvivw(i)
        resmed <- mr_mvmedian(i)
        resegger <- mvmr.eggerfx(bx, mrdat$b_out, bxse, mrdat$se_out)
        bind_rows(tibble(method = "ivw", exposure = trait1,
                         predictors = c(trait1, trait2),
                         b = resivw@Estimate, se = resivw@StdError,
                         ci_lo = resivw@CILower, ci_up = resivw@CIUpper,
                         p = resivw@Pvalue, nsnp = resivw@SNPs,
                         Q = resivw@Heter.Stat[1],
                         Q_df = resivw@SNPs - length(resivw@Exposure),
                         Q_pval = resivw@Heter.Stat[2]),
                  tibble(method = "wmedian", exposure = trait1,
                         predictors = c(trait1, trait2),
                         b = resmed@Estimate, se = resmed@StdError,
                         ci_lo = resmed@CILower, ci_up = resmed@CIUpper,
                         p = resmed@Pvalue, nsnp = resmed@SNPs),
                  mutate(resegger, exposure = trait1,
                         predictors = c("intercept", trait1, trait2))
                  ) %>%
            mutate(test = paste(trait1, "adj"))
    }) %>%
    bind_rows

rio::export(mvmr_res, "../files/mvmr_res.tsv")

mvmr_res %>%
    rename_with(toupper) %>%
    rename(CI_LOWER = CI_LO, CI_UPPER = CI_UP,
           BETA = B, MODEL = TEST) %>%
    mutate(METHOD = gsub("Ivw", "IVW", str_to_sentence(METHOD)),
           TRAIT = str_replace_all(
               EXPOSURE,
               c("sbp" = "SBP", "dbp" = "DBP",
                 "di" = "Insulin disposition index",
                 "intercept" = "Intercept", "whr" = "WHR",
                 "isi" = "Insulin sensitivity index")),
           MODEL = paste0(toupper(MODEL), "USTED"),
           across(c(BETA, SE, CI_LOWER, CI_UPPER), exp)) %>%
    relocate(MODEL) %>%
    rio::export("../files/mvmr_res.xlsx")
