## MR from exposure to T2D

## Libraries
library(tidyverse)

## Instruments
mr_ivs <- rio::import("../files/mr_ivs.tsv")

## Meta analysis of fat mass in legs and arms for each SNP
fm_df <- mr_ivs %>%
    filter(grepl("fatmass", trait_exp)) %>%
    mutate(trait_exp = gsub("_.+", "", trait_exp))

fm_meta <- fm_df %>%
    group_by(trait_exp, rsid) %>%
    summarise(beta.exposure = mean(beta.exposure),
              se.exposure = sqrt(sum(se.exposure^2) + (2 * 0.99 * prod(se.exposure))) / 2,
              eaf.exposure = mean(eaf.exposure), n.exposure = mean(n.exposure),
              .groups = "drop") %>%
    mutate(pval.exposure = 2 * pnorm(-abs(beta.exposure / se.exposure)))

fm <- fm_df %>%
    select(-names(fm_meta)[-c(1,2)]) %>%
    unique %>%
    inner_join(fm_meta)

mr_ivs <- mr_ivs %>%
    filter(!grepl("fatmass", trait_exp)) %>%
    bind_rows(fm)

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
    ## GWAS significant instruments for traits
    filter(rsid %in% snp_in_ref$V2, pval.exposure < 5e-8) %>%
    group_by(trait_exp) %>%
    group_modify(~{
        to_clump <- transmute(.x, rsid, pval = pval.exposure)
        clump_dat <- ieugwasr::ld_clump(to_clump, clump_kb = 500, clump_r2 = 0.01,
                                        plink_bin = genetics.binaRies::get_plink_binary(),
                                        bfile = "../files/mrtot2d_ref/ref")
        filter(.x, rsid %in% clump_dat$rsid)
    }) %>%
    ungroup %>%
    transmute(SNP = rsid,
              id.exposure = trait_exp, exposure = trait_exp,
              beta.exposure, se.exposure, pval.exposure, eaf.exposure,
              id.outcome = "t2d", outcome = "t2d",
              beta.outcome = beta.t2d, se.outcome = se.t2d,
              pval.outcome = p.t2d, mr_keep = TRUE)

rio::export(pheno_exp, "../files/pheno_exp.tsv")
