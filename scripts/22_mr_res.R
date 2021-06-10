## Preparing data for MR

library(tidyverse)
library(ieugwasr)

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
    "#Retic",        FALSE,           "ukb-d-30250_irnt", NA,
    "Platelet volume", FALSE, "ebi-a-GCST004599", NA,
    "CRP", TRUE, "ukb-d-30710_irnt", "ENSG00000132693" 
) %>%
    separate_rows(ensembl_id)

## Trait metadata (Downloaded in previous step)
traitstouse <- rio::import("~/dva/files/traitstouse.tsv")

exposures <- left_join(exposures, traitstouse) %>%
    select(-c(ncase, ncontrol, n_min, unit))

## Completing metadata
availg <- gwasinfo() %>%
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

exposures <- exposures %>% left_join(cisreg) %>% select(-trait)

## Extracting GWAS significant hits for each trait
## Note1: The MRC IEU database has aligned all effects to the forward strand
## Note2: Sample size and EAF have missing values
sigins <- tophits(unique(exposures$id), clump = 0)

rio::export(sigins, "~/dva/files/sigins.tsv")

sigins <- rio::import("~/dva/files/sigins.tsv")

## SNPs to be used as instruments
snp_ins <- sigins %>%
    ## Removing INDELs
    filter(nchar(ea) == 1,
           nchar(nea) == 1) %>%
    ## Setting key to join with outcome data
    mutate(chrom_pos = paste(chr, position, sep = ":")) %>%
    ## Removing possible multiallelic SNPs
    group_by(id, chrom_pos) %>%
    filter(n() == 1) %>%
    ungroup %>%
    ## Adding exposure metadata
    inner_join(exposures) %>%
    ## For traits requiring cis MR:
    mutate(mr_keep = ifelse(requireCisMR,
                            "&"(chr == chrom_gene,
                                "&"(position > minpos,
                                    position < maxpos)),
                            TRUE)) %>%
    filter(mr_keep) %>%
    ## Completing sample sizes
    mutate(n = coalesce(n, sample_size)) %>%
    select(-c(mr_keep, sample_size))

snp_ins <- rename_with(snp_ins, ~ paste(.x, "exposure", sep = "."),
                       c(ea, nea, eaf, beta, se, p, n, id))

rio::export(snp_ins, "~/dva/files/snp_ins.tsv")

snp_ins <- rio::import("~/dva/files/snp_ins.tsv")

## Significant variants for BMI
bmi_sig <- read_table2(pipe("awk 'NR==1||$9<5e-8' ~/dva/files/bmi.txt"))

bmi_sig <- select(bmi_sig, -c(CHR, POS))

names(bmi_sig) <- c("rsid",
                    paste(c("ea", "nea", "eaf", "beta", "se", "p", "n"), "bmi", sep = "."))

bmi_sig <- bmi_sig %>%
    mutate(ea = ifelse(beta.bmi > 0, ea.bmi, nea.bmi),
           nea = ifelse(beta.bmi > 0, nea.bmi, ea.bmi),
           beta.bmi = abs(beta.bmi)) %>%
    select(-c(ea.bmi, nea.bmi))

traitbmi_ins <- inner_join(snp_ins, bmi_sig)

traitbmi_ins <- traitbmi_ins %>%
    mutate(harmon = case_when("&"(ea == ea.exposure,
                                  nea == nea.exposure) ~ 1,
                              "&"(nea == ea.exposure,
                                  ea == nea.exposure) ~ -1,
                              TRUE ~ 0),
           beta.exposure = beta.exposure * harmon,
           eaf.exposure = ifelse(harmon == 1, eaf.exposure, 1 - eaf.exposure)) %>%
    filter(harmon != 0) %>%
    select(-c(harmon, ea.exposure, nea.exposure))

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
                                         paste(c("ea", "nea", "beta", "se",
                                                 "pval", "n"),
                                               "outcome", sep = "."))) %>%
    filter(nchar(ea.outcome) == 1,
           nchar(nea.outcome) == 1) %>%
    group_by(chrom_pos) %>%
    filter(n() == 1) %>%
    ungroup %>%    
    mutate(id.outcome = "t2d", outcome = "t2d")

## Harmonise data
harmondat_all <- inner_join(traitbmi_ins, outcome_t2d)

harmondat_all <- harmondat_all %>%
    mutate(harmon = case_when("&"(ea == ea.outcome,
                                  nea == nea.outcome) ~ 1,
                              "&"(ea == nea.outcome,
                                  nea == ea.outcome) ~ -1,
                              TRUE ~ 0),
           beta.outcome = beta.outcome * harmon) %>%
    filter(harmon != 0) %>%
    select(-c(harmon, ea.outcome, nea.outcome))

## Making reference panel to accelerate clumping
write(unique(harmondat_all$rsid), "~/dva/files/harmonsnps_all.txt", ncolumns = 1)
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
                      select = c(2,5)) %>% rename(rsid = SNP)

harmondat_all <- inner_join(harmondat_all, maf1kg)

## Calculating standardized effect sizes for the exposures
harmondat_all <- harmondat_all %>%
    mutate(z.exposure = beta.exposure / se.exposure,
           se.exposure = 1/sqrt(2 * MAF * (1 - MAF) * (n.exposure + (z.exposure^2))),
           beta.exposure = z.exposure * se.exposure) %>%
    select(-z.exposure)

rio::export(harmondat_all, "~/dva/files/harmondat_all.tsv")

## Clumping to obtain independent SNPs
mrdat_all <- harmondat_all %>%
    group_by(id.exposure) %>%
    group_modify(~{
        dat <- transmute(.x, rsid, pval = p.exposure)
        clump_dat <- ieugwasr::ld_clump(dat, clump_kb = 500, clump_r2 = 0.01,
                                        plink_bin = genetics.binaRies::get_plink_binary(),
                                        bfile = "/ludc/Home/daniel_c/dva/files/harmonsnps_all_ref/ref")
        filter(.x, rsid %in% clump_dat$rsid)
    }) %>%
    ungroup

rio::export(mrdat_all, "~/dva/files/mrdat_all.tsv")

mrdat_all <- rio::import("~/dva/files/mrdat_all.tsv")

## MR
source(list.files(pattern = "mrfx.R$"))

mr_res <- mrdat_all %>%
    group_by(trait_short) %>%
    filter(n() > 3) %>%
    group_modify(~mr_analysis(.x$beta.exposure, .x$beta.outcome, .x$se.exposure, .x$se.outcome))

rio::export(mr_res, "~/dva/files/mr_res.tsv")

## Traits with less than 3 variants
mrdat_all %>%
    group_by(trait_short) %>%
    filter(n() <= 3) %>%
    data.frame
