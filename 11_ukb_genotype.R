## UKBB genotype

library(tidyverse)

## BMI-T2D SNPs
mix <- rio::import("../files/mix.txt")

## Joining
grs_snps <- select(mix, rsid, chr)

## Creating a SNP file for PLINK per chromosome
if(!dir.exists("../files/grs_snps"))
    dir.create("../files/grs_snps")

system("rm ../files/grs_snps/*")

grs_snps %>%
    group_by(chr) %>%
    group_walk(~write(.x$rsid,
                      paste0("../files/grs_snps/rsids_chr",
                             unique(.y$chr), ".txt")))

### Querying UKBB

if(!dir.exists("../files/ukb_genotype"))
    dir.create("../files/ukb_genotype")

system("rm ../files/ukb_genotype/*")

chromosomes <- paste(unique(grs_snps$chr), collapse = " ")

system(
    paste0(
        "myarray=(", chromosomes,"); ",
        "for chr in \"${myarray[@]}\"; do ",
        "/ludc/Tools/Software/Plink/v2.00a1LM/plink2 ",
        "--bgen /ludc/Raw_Data_Archive/UKBB/imp/ukb_imp_chr${chr}_v3.bgen ",
        "--sample ../files/ukb_57232/ukb57232_imp_chr1_v3_s487296.sample ",
        "--extract ../files/grs_snps/rsids_chr${chr}.txt ",
        "--make-pgen ",
        "--out ../files/ukb_genotype/geno_chr$chr; done"
    )
)



