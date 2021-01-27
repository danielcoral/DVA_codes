## MR from exposure to T2D

## Libraries
library(tidyverse)

## Instruments
mr_ivs <- rio::import("../files/mr_ivs.tsv") %>%
    { ## Meta analysis of fat mass in legs and arms for each SNP
        df <- .
        fm_df <- filter(df, grepl("fatmass", trait_exp)) %>%
            mutate(trait_exp = gsub("_.+", "", trait_exp)) %>%
            group_by(trait_exp, rsid)
        fm_meta <- fm_df %>%
            summarise(b = mean(beta.exposure),
                      ## Assuming individuals in UKB are mostly symmetrical
                      se = sqrt(sum(se.exposure^2) + (2 * 0.99 * prod(se.exposure))) / 2,
                      z = b / se, p = 2 * pnorm(z, lower.tail = F),
                      eaf = mean(eaf.exposure), n = mean(n.exposure))
        fm_df %>%
            select(-c(beta.exposure, se.exposure, pval.exposure, eaf.exposure, n.exposure)) %>%
            unique %>%
            inner_join(rename(fm_meta, beta.exposure = b, se.exposure = se,
                              pval.exposure = p, eaf.exposure = eaf, n.exposure = n),
                       by = c("trait_exp", "rsid")) %>%
            bind_rows(filter(df, !grepl("fatmass", trait_exp))) }

## Making reference panel to accelerate clumping
write(unique(mr_ivs$rsid), "../files/mrtot2d_snps.txt", ncolumns = 1)
if(!dir.exists("../files/mrtot2d_ref/"))
    dir.create("../files/mrtot2d_ref/")
system("rm ../files/mrtot2d_ref/*")
system(paste0(genetics.binaRies::get_plink_binary(),
              " --bfile ../files/1kg_ref/EUR --extract ../files/mrtot2d_snps.txt ",
              "--make-bed --out ../files/mrtot2d_ref/ref"))

## Available SNPs
snp_in_ref <- rio::import("../files/mrtot2d_ref/ref.bim", format = "\t")

## Selecting instruments and parsing data
pheno_exp <- mr_ivs %>%
    ## GWAS significant instruments
    filter(rsid %in% snp_in_ref$V2, pval.exposure < 5e-8) %>%
    group_by(trait_exp, disc) %>%
    group_modify(~{
        to_clump <- transmute(.x, rsid, pval = pval.exposure)
        clump_dat <- ieugwasr::ld_clump(to_clump, clump_kb = 500, clump_r2 = 0.01,
                                        plink_bin = genetics.binaRies::get_plink_binary(),
                                        bfile = "../files/mrtot2d_ref/ref")
        filter(.x, rsid %in% clump_dat$rsid)
    }) %>%
    ungroup %>%
    transmute(SNP = rsid, disc,
              id.exposure = trait_exp, exposure = trait_exp, beta.exposure, se.exposure,
              id.outcome = "t2d", outcome = "t2d", beta.outcome = beta.t2d, se.outcome = se.t2d,
              mr_keep = T)

rio::export(pheno_exp, "../files/pheno_exp.tsv")

pheno_exp <- rio::import("../files/pheno_exp.tsv")

## MR analysis
## For this we will use three methods
## 1. IVW which handles the situation of valid IVs and balanced pleiotropy
## 2. MR-Egger which test if there is balanced pleiotropy
## (and if the dose-response relationship of IVW is correct)
## 3. Contamination-mixture which handles better directional (unbalanced) pleiotropy
## (pleiotropy through a confounder between the exposure and outcome)

## 1. Any evidence of an effect of the exposure on the outcome:
ivw_res <- pheno_exp %>%
    group_by(disc) %>%
    group_modify(~TwoSampleMR::mr(.x, method_list = "mr_ivw")) %>%
    mutate(method = "ivw", ci_lo = b - qnorm(1-0.05/2)*se, ci_up = b + qnorm(1-0.05/2)*se) %>%
    select(-c(id.exposure, id.outcome, outcome, se))

## Is the mean pleoiotropic effect = 0?
mr_pleio <- pheno_exp %>%
    group_by(disc) %>%
    group_modify(~bind_rows(TwoSampleMR::mr(.x, method_list = "mr_egger_regression") %>%
                            mutate(method = "egger"),
                            TwoSampleMR::mr_pleiotropy_test(.x) %>%
                            rename(b = egger_intercept) %>%
                            mutate(method = "intercept"))) %>%
    mutate(ci_lo = b - qnorm(1-0.05/2)*se, ci_up = b + qnorm(1-0.05/2)*se) %>%
    select(-c(id.exposure, id.outcome, outcome, se))

## Do the results hold in unbalanced pleoitropy?
library(MendelianRandomization)

conmix_res <- pheno_exp %>%
    group_by(exposure, disc) %>%
    group_modify(~{
        lapply(seq(1,2,.1),
               function(psi){
                   res <- mr_input(bx = .x$beta.exposure, bxse = .x$se.exposure,
                                   by = .x$beta.outcome, byse = .x$se.outcome) %>%
                       mr_conmix(psi = psi * sd(.x$beta.outcome/.x$beta.exposure))
                   data.frame(b = res@Estimate, ci_lo = res@CILower, ci_up = res@CIUpper,
                              pval = res@Pvalue, nsnp = res@SNPs, psi = psi,
                              psitrue = res@Psi, method = "conmix")
               }) %>%
            bind_rows()
    })

mr_res <- bind_rows(ivw_res, mr_pleio, conmix_res)

rio::export(mr_res, "../files/mr_res.tsv")
