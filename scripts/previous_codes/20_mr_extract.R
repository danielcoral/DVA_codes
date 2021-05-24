## MR data extraction

library(tidyverse)

## BMI - T2D instruments
sigbmi_t2d <- rio::import("../files/sigbmi_t2d.RData")

## Data from MRC IEU database
ids_ieu <- c(
    ## The betas of all these variables are in SD units
    armfatmass_l = "ukb-b-8338", armfatmass_r = "ukb-b-6704",
    legfatmass_l = "ukb-b-7212", legfatmass_r = "ukb-b-18096",
    hdl = "ukb-d-30760_irnt", ldl = "ukb-d-30780_irnt",
    tg = "ukb-d-30870_irnt"
)

## Extracting instrument from MRC IEU data
gwasvcf::set_bcftools()
ieu_dat <- lapply(ids_ieu,
                  function(x){
                      message(paste("Extracting variants from", x, "-",
                                    which(ids_ieu == x), "of", length(ids_ieu), "..."),
                              appendLF = FALSE)
                      res <- gwasvcf::query_gwas(paste0("../files/mr_ss/", x, ".vcf.gz"),
                                                 chrompos = sigbmi_t2d$chr_pos) %>%
                          ## Converting data from gwasvcf
                          gwasglue::gwasvcf_to_TwoSampleMR(type="exposure") %>%
                          ## Parsing column names
                          rename(rsid = SNP, chr = chr.exposure, pos = pos.exposure,
                                 ea.exposure = effect_allele.exposure,
                                 nea.exposure = other_allele.exposure,
                                 n.exposure = samplesize.exposure) %>%
                          mutate(chr = as.numeric(as.character(chr))) %>%
                          select(-paste0(c("ncase.", "ncontrol.", "mr_keep.",
                                           "pval_origin.", "id.", ""),
                                         "exposure")) %>%
                          ## Adding T2D information
                          inner_join(sigbmi_t2d, by = c("rsid", "chr", "pos")) %>%
                          ## Only one association for each original SNP
                          group_by(orig_rsid) %>%
                          arrange(is_proxy, r2, pval.exposure) %>%
                          slice(1) %>%
                          ungroup() %>%
                          ## Aligning to BMI increasing allele
                          mutate(harmon = case_when(
                                     ea.exposure == ea & nea.exposure == nea ~ 1,
                                     ea.exposure == nea & nea.exposure == ea ~ -1,
                                     T ~ 0),
                                 beta.exposure = harmon * beta.exposure) %>%
                          filter(harmon != 0) %>%
                          select(-c(harmon, chr_pos))
                      message("done")
                      return(res)
                  })

ieu_dfs <- bind_rows(ieu_dat, .id = "trait_exp")

## Blood pressure
bp_dfs <- lapply(
    c("sbp", "dbp"),
    function(x){
        vroom::vroom(paste0("../files/mr_ss/", x, ".txt.gz")) %>%
            `names<-`(c("chr_pos", "ea.exposure", "nea.exposure", "eaf.exposure",
                        "beta.exposure", "se.exposure", "pval.exposure",
                        "n.exposure", "neff")) %>%
            ## Removing indels
            filter(!grepl("INDEL", chr_pos)) %>%
            mutate(chr_pos = gsub(":SNP", "", chr_pos), trait_exp = x,
                   ea.exposure = toupper(ea.exposure),
                   nea.exposure = toupper(nea.exposure),
                   neff = NULL) %>%
            inner_join(sigbmi_t2d, by = "chr_pos") %>%
            separate(chr_pos, c("chr", "pos"), sep = ":", convert = T) %>%
            group_by(orig_rsid) %>%
            arrange(is_proxy, r2, pval.exposure) %>%
            slice(1) %>%
            ungroup() %>%
            mutate(harmon = case_when(
                       ea.exposure == ea & nea.exposure == nea ~ 1,
                       ea.exposure == nea & nea.exposure == ea ~ -1,
                       T ~ 0),
                   beta.exposure = harmon * beta.exposure) %>%
            filter(harmon != 0) %>%
            ## Finding standardized effect size (Zhu et al. 2018)
            mutate(zscore = beta.exposure / se.exposure,
                   beta.exposure = zscore /
                       sqrt(2 * eaf.exposure * (1-eaf.exposure) * (n.exposure + zscore^2)),
                   se.exposure = 1 /
                       sqrt(2 * eaf.exposure * (1-eaf.exposure) * (n.exposure + zscore^2))) %>%
            select(-c(harmon, zscore))
    }
)

bp_dfs <- bind_rows(bp_dfs)

## WHR
whr <- vroom::vroom("../files/mr_ss/whr.txt.gz") %>%
    `names<-`(c("chr", "pos", "rsid", "ea.exposure", "nea.exposure", "eaf.exposure",
                "beta.exposure", "se.exposure", "pval.exposure", "n.exposure", "INFO")) %>%
    ## Filter for INFO score and MAF > 0.01
    filter(INFO > 0.8, eaf.exposure > 0.01, eaf.exposure < 0.99) %>%
    mutate(rsid = gsub(":.+", "", rsid)) %>%
    ## No need to filter INDELS, they are not included
    ## Adding T2D information
    inner_join(sigbmi_t2d, by = c("rsid", "chr", "pos")) %>%
    ## Only one association for each original SNP
    group_by(orig_rsid) %>%
    arrange(is_proxy, r2, pval.exposure) %>%
    slice(1) %>%
    ungroup() %>%
    ## Aligning to BMI increasing allele
    mutate(harmon = case_when(
               ea.exposure == ea & nea.exposure == nea ~ 1,
               ea.exposure == nea & nea.exposure == ea ~ -1,
               T ~ 0),
           beta.exposure = harmon * beta.exposure,
           trait_exp = "whr") %>%
    filter(harmon != 0) %>%
    select(-c(harmon, chr_pos, INFO))

mr_ivs <- bind_rows(ieu_dfs, bp_dfs, whr)

rio::export(mr_ivs, "../files/mr_ivs.tsv")
