library(readr)
library(dplyr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(ggplot2)
bmi <- read_tsv("~/projects/DVA/Data/GWAS_sumstats/bmi.txt", 
                skip = 1,
                col_names = c("chrom", "pos", "rsid", 
                              paste(c("ea", "nea", "eaf", "beta", "se", "pval", "n"), "bmi", sep = "_")),
                col_types = "nncccnnnnn")
head(bmi)
t2d <- read_tsv("~/projects/DVA/Data/GWAS_sumstats/t2d.txt", 
                skip = 1,
                col_names = c("SNP", "chrom", "pos",
                              paste(c("ea", "nea", "eaf", "beta", "se", "pval", "n"), "t2d", sep = "_")),
                col_types = "cnnccnnnnn")
head(t2d)
bmi_t2d <- bmi %>%
    unite(SNP, chrom, pos, sep = ":") %>%
    inner_join(t2d, by = "SNP") %>%
    select(-SNP) %>%
    ## MAF concordance (Difference less than 20%)
    mutate(maf_bmi = ifelse(eaf_bmi < .5, eaf_bmi, 1 - eaf_bmi),
           maf_t2d = ifelse(eaf_t2d < .5, eaf_t2d, 1 - eaf_t2d)) %>%
    filter(abs(maf_bmi - maf_t2d) < .2) %>%
    ## Removing palindromic SNPs with MAF > 40%
    mutate(palind1_bmi = ea_bmi %in% c("A", "T") & nea_bmi %in% c("A", "T"),
           palind2_bmi = ea_bmi %in% c("C", "G") & nea_bmi %in% c("C", "G"),
           palind1_t2d = ea_t2d %in% c("A", "T") & nea_t2d %in% c("A", "T"),
           palind2_t2d = ea_t2d %in% c("C", "G") & nea_t2d %in% c("C", "G"),
           palind_any = palind1_bmi | palind2_bmi | palind1_t2d | palind2_t2d,
           palind_ambig = palind_any & (maf_bmi > .4 | maf_t2d > .4)) %>%
    filter(!palind_ambig) %>%
    select(-starts_with("palind"), -maf_bmi, -maf_t2d) %>%
    ## Aligning to the BMI increasing allele
    mutate(ea = ifelse(beta_bmi >= 0, ea_bmi, nea_bmi),
           nea = ifelse(beta_bmi >= 0, nea_bmi, ea_bmi),
           eaf_bmi = ifelse(beta_bmi >= 0, eaf_bmi, 1 - eaf_bmi),
           beta_bmi = abs(beta_bmi)) %>%
    select(-c(ea_bmi, nea_bmi)) %>%
    mutate(harmon = case_when(ea == ea_t2d & nea == nea_t2d ~ 1,
                              ea == nea_t2d & nea == ea_t2d ~ -1,
                              TRUE ~ 0),
           ## If the alleles don't match, flipping strand
           across(c(ea_t2d, nea_t2d),
                  ~ifelse(harmon == 0,
                          toupper(
                              stringr::str_replace_all(.x, c("A" = "t", "T" = "a", "C" = "g", "G" = "c"))
                              ),
                          .x)),
           ## Retesting with flipped alleles
           harmon = ifelse(harmon == 0,
                           case_when(harmon == 0 & ea == ea_t2d & nea == nea_t2d ~ 1,
                                     harmon == 0 & ea == nea_t2d & nea == ea_t2d ~ -1,
                                     TRUE ~ 0),
                           harmon),
           ## Aligning
           beta_t2d = beta_t2d * harmon,
           eaf_t2d = ifelse(harmon == 1, eaf_t2d, 1 - eaf_t2d)) %>%
    ## Removing SNPs with allele mismatch
    filter(harmon != 0) %>%
    select(-c(harmon, ea_t2d, nea_t2d)) %>%
    select(chrom, pos, rsid, ea, nea, everything()) %>%
    arrange(chrom, pos)
head(bmi_t2d)
nrow(bmi_t2d)
write_tsv(bmi_t2d, "../data/bmi_t2d.tsv")
fn <- tempfile()
bmi_t2d %>%
    filter(pval_bmi < 5e-8, pval_t2d < 5e-8) %>%
    transmute(SNP = rsid, P = pval_bmi) %>%
    write_tsv(fn)
plink_cmd <- paste0("/ludc/Tools/Software/Plink/v1.90b5.2/plink",
                    " --bfile ", "~/projects/DVA/Data/ReferenceData/1kg_ref/EUR",
                    " --clump ", fn, 
                    " --clump-p1 ", 5e-8, 
                    " --clump-r2 ", 0.01,
                    " --clump-kb ", 500, 
                    " --out ", fn)
system(plink_cmd)
clumped_snps <- paste("awk '{print $3}'", paste(fn, "clumped", sep = ".")) %>%
    pipe %>%
    readLines
unlink(paste0(fn, "*"))
head(clumped_snps)
bmi_t2d_clumped <- filter(bmi_t2d, rsid %in% clumped_snps)
head(bmi_t2d_clumped)
nrow(bmi_t2d_clumped)
gencode <- rtracklayer::readGFF("~/projects/DVA/Data/ReferenceData/gencode.v19.annotation.gtf")
head(gencode)
bmi_t2d_genprofiles <- gencode %>%
    filter(type == "gene") %>%
    transmute(chrom = as.numeric(gsub("chr", "", seqid)),
              start, end, gene_name) %>%
    drop_na %>%
    inner_join(bmi_t2d_clumped, by = "chrom") %>%
    mutate(distance = ifelse(pos >= start & pos <= end, 0, pmin(abs(pos - start), abs(pos - end)))) %>%
    group_by(rsid) %>%
    slice_min(distance, with_ties = FALSE) %>%
    mutate(profile = ifelse(beta_t2d > 0, "Concordant", "Discordant"), nearest_gene = gene_name) %>%
    select(-c(start, end, gene_name, distance))
head(bmi_t2d_genprofiles)
write_tsv(bmi_t2d_genprofiles, "../data/mix.tsv")
