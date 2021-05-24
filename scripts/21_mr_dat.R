## Preparing data for MR

library(tidyverse)
library(TwoSampleMR)

## Relevant exposures identified in the previous analyses
exposures <- tribble(
    ~ trait_short,   ~ requireCisMR,  ~ id,               ~ ensembl_id,
    "SHBG",          TRUE,            "ukb-d-30830_irnt", "ENSG00000129214",
    "Urate",         FALSE,           "ukb-d-30880_irnt", NA,
    "GGT",           TRUE,            "ukb-d-30730_irnt", "ENSG00000149435;ENSG00000100121;ENSG00000274252;ENSG00000230712;ENSG00000276160;ENSG00000100031;ENSG00000133475;ENSG00000197421;ENSG00000280208;ENSG00000099998;ENSG00000167741;ENSG00000131067;ENSG00000236969",
    "WHR",           FALSE,           "ieu-a-73",         NA,         
    "AdipoQ",        TRUE,            "ieu-a-1",          "ENSG00000181092",
    "ALT",           TRUE,            "ukb-d-30620_irnt", "ENSG00000167701;ENSG00000166123",
    "AST",           TRUE,            "ukb-d-30650_irnt", "ENSG00000120053;ENSG00000125166",
    "#Leuco",        FALSE,           "ukb-d-30000_irnt", NA,
    "SBP",           FALSE,           "ieu-b-38",         NA,
    "DBP",           FALSE,           "ieu-b-39",         NA,
    "HDL",           FALSE,           "ukb-d-30760_irnt", NA,
    "ApoA1",         TRUE,            "ieu-b-107",        "ENSG00000118137",
    "LDL",           FALSE,           "ukb-d-30780_irnt", NA,
    "TG",            FALSE,           "ukb-d-30870_irnt", NA,
    "FVC",           FALSE,           "ukb-b-7953",       NA,
    "HeelBMD",       FALSE,           "ukb-b-20124",      NA,
    "PTH",           TRUE,            "prot-a-2431",      "ENSG00000152266",
    "Phosphate",      FALSE,           "ukb-d-30810_irnt", NA,
    "TotalLeanMass", FALSE,           "ukb-b-13354",      NA,
    "ArmFatMassR",   FALSE,           "ukb-b-6704",       NA,
    "LegFatMassR",   FALSE,           "ukb-b-18096",      NA,
    "#Retic",        FALSE,           "ukb-d-30250_irnt", NA
) %>%
    separate_rows(ensembl_id)

## Trait metadata (Downloaded in previous step)
traitstouse <- rio::import("~/dva/files/traitstouse.tsv")

exposures <- left_join(exposures, traitstouse) %>%
    select(-c(ncase, ncontrol, n_min, unit))

## Completing metadata
availg <- available_outcomes() %>%
    select(all_of(intersect(names(exposures), names(.))))

missingdat <- exposures %>%
    filter(is.na(category)) %>%
    select(trait_short, requireCisMR, id, ensembl_id)

toreplace <- inner_join(missingdat, availg) %>%
    mutate(category = "Continuous", sex = "Both")

exposures <- exposures %>%
    filter(!id %in% toreplace$id) %>%
    bind_rows(toreplace)

## Identifying cis-regions using Ensembl BiomaRt
ensembl_genes <- biomaRt::useEnsembl(biomart = "ensembl",
                                     dataset = "hsapiens_gene_ensembl",
                                     GRCh = 37)

ensgids <- exposures %>% filter(!is.na(ensembl_id)) %>% pull(ensembl_id)

cisreg <- biomaRt::getBM(attributes = c("ensembl_gene_id", "chromosome_name",
                                        "start_position", "end_position",
                                        "transcription_start_site"),
                         filters = "ensembl_gene_id",
                         values = unique(na.exclude(exposures$ensembl_id)),
                         mart = ensembl_genes) %>%
    group_by(ensembl_gene_id) %>%
    summarise(chrom_gene = unique(chromosome_name),
              minpos = min(min(start_position),
                           min(end_position),
                           min(transcription_start_site)) - 1000000,
              maxpos = max(max(start_position),
                           max(end_position),
                           max(transcription_start_site)) + 1000000) %>%
    rename(ensembl_id = ensembl_gene_id)

exposures <- exposures %>% left_join(cisreg)

## Extracting GWAS significant hits for each trait
## Note: The MRC IEU database has aligned all effects to the forward strand
sigins <- map_dfr(unique(exposures$id), ~ extract_instruments(.x, clump = FALSE))

rio::export(sigins, "~/dva/files/sigins.tsv")

sigins <- rio::import("~/dva/files/sigins.tsv")

## SNPs to be used as instruments
snp_ins <- sigins %>%
    ## Removing INDELs
    filter(nchar(effect_allele.exposure) == 1,
           nchar(other_allele.exposure) == 1) %>%
    ## Setting key to join with outcome data
    mutate(chrom_pos = paste(chr.exposure, pos.exposure, sep = ":")) %>%
    ## Removing possible multiallelic SNPs
    group_by(id.exposure, chrom_pos) %>%
    filter(n() == 1) %>%
    ungroup %>%
    ## Adding exposure metadata
    inner_join(exposures, by = c("id.exposure" = "id")) %>%
    ## For traits requiring cis MR:
    mutate(mr_keep = ifelse(
               requireCisMR,
               "&"(chr.exposure == chrom_gene,
                   "&"(pos.exposure > minpos,
                       pos.exposure < maxpos)),
               mr_keep.exposure
           )) %>%
    filter(mr_keep) %>%
    ## Completing sample sizes
    mutate(samplesize.exposure = coalesce(samplesize.exposure,
                                          sample_size)) %>%
    select(-c(mr_keep.exposure, sample_size))

rio::export(snp_ins, "~/dva/files/snp_ins.tsv")

## Outcome: T2D from DIAGRAM
## Note: These results are aligned to the forward strand

## Subsetting SNPs needed
snp_ins %>%
    select(chrom_pos) %>%
    unique %>%
    rio::export("~/dva/files/snp_ins_chrompos.tsv")

system(paste("grep -Fwf",
             "~/dva/files/snp_ins_chrompos.tsv",
             "~/dva/files/t2d_diagram.txt",
             "> ~/dva/files/outcome_t2d.tsv"))

outcome_t2d <- rio::import("~/dva/files/outcome_t2d.tsv",
                           col.names = c("chrom_pos",
                                         paste(c("effect_allele",
                                                 "other_allele", "beta", "se",
                                                 "pval", "samplesize"),
                                               "outcome", sep = "."))) %>%
    filter(nchar(effect_allele.outcome) == 1,
           nchar(other_allele.outcome) == 1) %>%
    group_by(chrom_pos) %>%
    filter(n() == 1) %>%
    ungroup %>%    
    mutate(id.outcome = "t2d", outcome = "t2d")

## Harmonise data
harmondat_all <- inner_join(snp_ins, outcome_t2d) %>%
    ## Aligning to the allele increasing exposure
    transmute(
        chrom = chr.exposure, pos = pos.exposure, SNP,
        effect_allele = ifelse(beta.exposure > 0,
                               effect_allele.exposure,
                               other_allele.exposure),
        other_allele = ifelse(beta.exposure > 0,
                              other_allele.exposure,
                              effect_allele.exposure),
        id.exposure, exposure = trait_short,
        beta.exposure = abs(beta.exposure), se.exposure,
        pval.exposure, samplesize.exposure,
        id.outcome = "t2d_diagram", outcome = "t2d",
        harmon = case_when(
            "&"(effect_allele == effect_allele.outcome,
                other_allele == other_allele.outcome) ~ 1,
            "&"(effect_allele == other_allele.outcome,
                other_allele == effect_allele.outcome) ~ -1,
            TRUE ~ 0
        ),
        beta.outcome = beta.outcome * harmon, se.outcome,
        pval.outcome, samplesize.outcome,
        mr_keep = TRUE,
        year, consortium, author, pmid, population,
        nsnp_total = nsnp
    ) %>%
    filter(harmon != 0) %>%
    select(-harmon)

## Making reference panel to accelerate clumping
write(unique(harmondat_all$SNP), "~/dva/files/harmonsnps_all.txt", ncolumns = 1)
if(!dir.exists("~/dva/files/harmonsnps_all_ref/"))
    dir.create("~/dva/files/harmonsnps_all_ref/")
system("rm ~/dva/files/harmonsnps_all_ref/*")
system(paste0(genetics.binaRies::get_plink_binary(),
              " --bfile ~/dva/files/1kg_ref/EUR",
              " --extract ~/dva/files/harmonsnps_all.txt",
              " --make-bed --out ~/dva/files/harmonsnps_all_ref/ref"))

## Extracting minor allele frequencies
system(paste(genetics.binaRies::get_plink_binary(),
             "--bfile /ludc/Home/daniel_c/dva/files/harmonsnps_all_ref/ref",
             "--freq",
             "--out /ludc/Home/daniel_c/dva/files/mrsnpfreq"))

maf1kg <- rio::import("/ludc/Home/daniel_c/dva/files/mrsnpfreq.frq", format = "\t",
                      select = c(2,5))

harmondat_all <- inner_join(harmondat_all, maf1kg)

## Calculating standardized effect sizes for the exposures
harmondat_all <- harmondat_all %>%
    mutate(z.exposure = beta.exposure / se.exposure,
           se.exposure = 1/sqrt(2 * MAF * (1 - MAF) * (samplesize.exposure + (z.exposure^2))),
           beta.exposure = z.exposure * se.exposure) %>%
    select(-z.exposure)

rio::export(harmondat_all, "~/dva/files/harmondat_all.tsv")

## Clumping to obtain independent SNPs
mrdat_all <- harmondat_all %>%
    group_by(id.exposure) %>%
    group_modify(~{
        dat <- transmute(.x, rsid = SNP, pval = pval.exposure)
        clump_dat <- ieugwasr::ld_clump(dat, clump_kb = 500, clump_r2 = 0.01,
                                        plink_bin = genetics.binaRies::get_plink_binary(),
                                        bfile = "/ludc/Home/daniel_c/dva/files/harmonsnps_all_ref/ref")
        filter(.x, SNP %in% clump_dat$rsid)
    }) %>%
    ungroup

rio::export(mrdat_all, "~/dva/files/mrdat_all.tsv")

## Obtaining BMI effects
system(paste("grep -Fwf",
             "~/dva/files/harmonsnps_all.txt",
             "~/dva/files/bmi.txt",
             "> ~/dva/files/bmi_ins.tsv"))

bmi_ins <- rio::import("~/dva/files/bmi_ins.tsv",
                       col.names = c("chrom", "pos", "SNP",
                                     paste(c("effect_allele", "other_allele",
                                             "eaf", "beta", "se",
                                             "pval", "samplesize"),
                                           "bmi", sep = "."))) %>%
    filter(nchar(effect_allele.bmi) == 1,
           nchar(other_allele.bmi) == 1) %>%
    group_by(SNP) %>%
    filter(n() == 1) %>%
    ungroup %>%
    mutate(chrom = as.character(chrom))

harmondat_bmi <- inner_join(harmondat_all, bmi_ins)

rio::export(harmondat_bmi, "~/dva/files/harmondat_bmi.tsv")

## For MR restricted to BMI significant variants
mrdat_bmisig <-  harmondat_bmi %>%
    filter(pval.bmi < 5e-8) %>%
    ## Aligning to the BMI increasing alleles
    mutate(ea.bmi = ifelse(beta.bmi > 0, effect_allele.bmi,
                           other_allele.bmi),
           nea.bmi = ifelse(beta.bmi > 0, other_allele.bmi,
                            effect_allele.bmi),
           beta.bmi = abs(beta.bmi),
           harmon = case_when("&"(effect_allele == ea.bmi,
                                  other_allele == nea.bmi) ~ 1,
                              "&"(effect_allele == nea.bmi,
                                  other_allele == ea.bmi) ~ -1,
                              TRUE ~ 0),
           beta.exposure = beta.exposure * harmon,
           beta.outcome = beta.outcome * harmon,
           effect_allele = ea.bmi, other_allele = nea.bmi) %>%
    filter(harmon != 0) %>%
    select(-c(effect_allele.bmi, other_allele.bmi,
              eaf.bmi, ea.bmi, nea.bmi, harmon)) %>%
    group_by(id.exposure) %>%
    group_modify(~{
        dat <- transmute(.x, rsid = SNP, pval = pval.exposure)
        clump_dat <- ieugwasr::ld_clump(dat, clump_kb = 500, clump_r2 = 0.01,
                                        plink_bin = genetics.binaRies::get_plink_binary(),
                                        bfile = "/ludc/Home/daniel_c/dva/files/harmonsnps_all_ref/ref")
        filter(.x, SNP %in% clump_dat$rsid)
    }) %>%
    ungroup

rio::export(mrdat_bmisig, "~/dva/files/mrdat_bmisig.tsv")
