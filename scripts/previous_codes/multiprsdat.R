## Finding subtypes of obesity in UKB

library(tidyverse)
library(data.table)
library(ieugwasr)

## Level of significance for the 4 groups of traits
## BP, lipids and central adiposity
pthres <- 0.01

## This is because
pthres ^ 4

## SNPs significant for BMI
bmi <- fread("../files/bmi.txt")

bmisig <-  bmi[
   ,.(chr = CHR, hg19_pos = POS, rsid = SNP,
      ea = Tested_Allele, nea = Other_Allele,
      eaf = Freq_Tested_Allele_in_HRS,
      beta = BETA, se = SE, p = P, n = N,
      trait = "bmi")
][
    p < pthres & eaf > 0.01 & eaf < 0.99
][
  , `:=`(pal1 = ea %in% c("A", "T") & nea %in% c("A", "T"),
         pal2 = ea %in% c("C", "G") & nea %in% c("C", "G"),
         maf = ifelse(eaf > 0.5, 1 - eaf, eaf))
][
    !((pal1 | pal2) & maf > 0.3)
][
  , `:=`(ea2 = ifelse(beta > 0, ea, nea),
         nea2 = ifelse(beta > 0, nea, ea))
][
  , `:=`(ea = ea2, nea = nea2, beta = abs(beta))
][
  , -(pal1:nea2)
]

## WHR GIANT + UKB (not available in MRC IEU database)
whr <- fread("../files/mr_ss/whr.txt.gz")

whrsig <- whr[
    INFO > 0.8 & P < pthres
][
   , .(chr = CHR, hg19_pos = POS,
       rsid = gsub(":.*", "", SNP),
       ea = Tested_Allele, nea = Other_Allele,
       eaf = Freq_Tested_Allele, beta = BETA, se = SE, p = P, n = N,
       trait = "whr")
][
    rsid %in% bmisig$rsid
]

snps <- bmisig %>%
    filter(rsid %in% whrsig$rsid) %>%
    mutate(MarkerName = paste(chr, hg19_pos, "SNP", sep = ":"))

## BP from UKB + ICBP
sbp <- fread("../files/mr_ss/sbp.txt.gz")
dbp <- fread("../files/mr_ss/dbp.txt.gz")

sbp_f <- sbp[ P < pthres ][ snps[,.(MarkerName, rsid)], on = "MarkerName", nomatch = 0 ] 

dbp_f <- dbp[ P < pthres ][ sbp_f[,.(MarkerName, rsid)], on = "MarkerName", nomatch = 0 ] 

## Lipids
gwasvcf::set_bcftools()

hdl <- gwasvcf::query_gwas("../files/mr_ss/ukb-d-30760_irnt.vcf.gz",
                              chrompos = gsub(":SNP", "", dbp_f$MarkerName)) %>%
    gwasglue::gwasvcf_to_TwoSampleMR(type="exposure") %>%
    filter(pval.exposure < pthres)

tg <- gwasvcf::query_gwas("../files/mr_ss/ukb-d-30870_irnt.vcf.gz",
                          chrompos = paste(hdl$chr.exposure, hdl$pos.exposure, sep = ":")) %>%
    gwasglue::gwasvcf_to_TwoSampleMR(type="exposure") %>%
    filter(pval.exposure < pthres)

mprsdat <- bmisig %>%
    filter(rsid %in% tg$SNP) %>%
    bind_rows(filter(whrsig, rsid %in% tg$SNP)) %>%
    bind_rows(sbp_f %>%
              separate(MarkerName, c("chr", "hg19_pos"), extra = "drop") %>%
              transmute(chr = as.numeric(chr), hg19_pos = as.numeric(hg19_pos),
                        rsid, ea = toupper(Allele1), nea = toupper(Allele2), trait = "sbp",
                        eaf = Freq1, beta = Effect, se = StdErr, p = P, n = TotalSampleSize) %>%
              filter(rsid %in% tg$SNP)) %>%
    bind_rows(dbp_f %>%
              separate(MarkerName, c("chr", "hg19_pos"), extra = "drop") %>%
              transmute(chr = as.numeric(chr), hg19_pos = as.numeric(hg19_pos),
                        rsid, ea = toupper(Allele1), nea = toupper(Allele2), trait = "dbp",
                        eaf = Freq1, beta = Effect, se = StdErr, p = P, n = TotalSampleSize) %>%
              filter(rsid %in% tg$SNP)) %>%
    bind_rows(hdl %>%
              rename_with(~gsub("\\.exposure", "", .x)) %>%
              rename(hg19_pos = pos, ea = effect_allele, nea = other_allele,
                     p = pval, n = samplesize, rsid = SNP) %>%
              mutate(trait = "hdl",
                     across(c(chr, hg19_pos), as.numeric)) %>%
              select(-c(ncase, ncontrol, exposure, mr_keep, pval_origin, id)) %>%
              filter(rsid %in% tg$SNP)) %>%
    bind_rows(tg %>%
              rename_with(~gsub("\\.exposure", "", .x)) %>%
              rename(hg19_pos = pos, ea = effect_allele, nea = other_allele,
                     p = pval, n = samplesize, rsid = SNP) %>%
              mutate(trait = "tg",
                     across(c(chr, hg19_pos), as.numeric)) %>%
              select(-c(ncase, ncontrol, exposure, mr_keep, pval_origin, id))) %>%
    group_by(rsid) %>%
    mutate(pval = p[trait == "bmi"],
           ea_bmi = ea[trait == "bmi"],
           nea_bmi = nea[trait == "bmi"]) %>%
    ungroup %>%
    mutate(beta = case_when(ea == ea_bmi & nea == nea_bmi ~ beta,
                            ea == nea_bmi & nea == ea_bmi ~ -beta,
                            TRUE ~ NaN)) %>%
    transmute(rsid, chr, hg19_pos, ea_bmi, nea_bmi, pval, trait, beta) %>%
    pivot_wider(names_from = trait, values_from = beta) %>%
    ld_clump(clump_kb = 500, clump_r2 = 0.01,
             bfile = "../files/1kg_ref/EUR",
             plink_bin = genetics.binaRies::get_plink_binary()) %>%
    mutate(snpgroup = case_when(sbp > 0 & dbp > 0 & tg > 0 & hdl < 0 & whr > 0 ~ "1",
                                sbp > 0 & dbp > 0 & tg > 0 & hdl < 0 & whr < 0 ~ "2",
                                sbp > 0 & dbp > 0 & tg < 0 & hdl > 0 & whr > 0 ~ "3",
                                sbp > 0 & dbp > 0 & tg < 0 & hdl > 0 & whr < 0 ~ "4",
                                sbp < 0 & dbp < 0 & tg > 0 & hdl < 0 & whr > 0 ~ "5",
                                sbp < 0 & dbp < 0 & tg > 0 & hdl < 0 & whr < 0 ~ "6",
                                sbp < 0 & dbp < 0 & tg < 0 & hdl > 0 & whr > 0 ~ "7",
                                sbp < 0 & dbp < 0 & tg < 0 & hdl > 0 & whr < 0 ~ "8"))

## Adding T2D information

t2d_diagram <- fread("../files/t2d_diagram.txt")

snps_out <- mprsdat %>%
    transmute(rsid, chrpos = paste(chr, hg19_pos, sep = ":")) %>%
    setDT

t2d_diagram[,chrpos := `Chr:Position`]

t2d_f <- t2d_diagram[snps_out, on = "chrpos", nomatch = 0]

gres <- t2d_f %>%
    transmute(rsid, ea_t2d = Allele1, nea_t2d = Allele2,
              beta_t2d = Effect, se_t2d = StdErr) %>%
    inner_join(mprsdat) %>%
    mutate(beta_t2d = case_when(ea_bmi == ea_t2d & nea_bmi == nea_t2d ~ beta_t2d,
                                ea_bmi == ea_t2d & nea_bmi == nea_t2d ~ -beta_t2d,
                                TRUE ~ NaN)) %>%
    select(-c(ea_t2d, nea_t2d)) %>%
    group_by(snpgroup) %>%
    group_modify(~{
        ## Within each stratum, calculating the average logOR
        res <- meta::metagen(TE = beta_t2d, seTE = se_t2d,
                             data = .x, comb.fixed = FALSE,
                             comb.random = TRUE, method.tau = "PM",
                             hakn = FALSE, prediction = FALSE)
        with(res, data.frame(b = TE.random, se = seTE.random,
                             ci_lo = lower.random, ci_up = upper.random,
                             p = pval.random, Qv = Q, Qdf = df.Q, Qp = pval.Q,
                             Tau = tau))
    }) %>%
    filter(Qdf != 0, !is.na(snpgroup))

gres %>%
    ggplot(aes(snpgroup, b)) +
    geom_point() +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_up))

        
