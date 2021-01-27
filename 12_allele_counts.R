## Preprocessing allele counts

library(tidyverse)

## Origin of SNPs counted in UKBB:

## BMI-T2D SNPs
mix <- rio::import("../files/mix.txt")

## Obtaining allele counts
if(!dir.exists("../files/allele_counts"))
    dir.create("../files/allele_counts")

system("rm ../files/allele_counts/*")

chromosomes <- paste(unique(mix$chr), collapse = " ")

system(
    paste0(
        "myarray=(", chromosomes, "); ",
        "for chr in \"${myarray[@]}\"; do ",
        "/ludc/Tools/Software/Plink/v2.00a2LM/plink2 ",
        "--pfile ../files/ukb_genotype/geno_chr${chr} ",
        "--export A ",
        "--out ../files/allele_counts/counts_$chr; ",
        "done"
    )
)

## SNP count in UKBB
snp_count <- lapply(
    ## Looping over files
    list.files("../files/allele_counts/",
               pattern = "raw", full.names = T),
    function(x){
        ## Importing
        rio::import(x, format = "\t") %>%
            ## Dropping unused columns
            select(-c(SEX, FID, MAT, PAT, PHENOTYPE))
    }
) %>%
    ## Joining all columns using the identifier
    reduce(inner_join, by = "IID") %>%
    ## Calculating absolute allele counts
    mutate_at(
        vars(starts_with("rs")), 
        function(x){
            case_when(
                x < 0.5 ~ 0,
                x >= 0.5 & x < 1.5 ~ 1,
                x >= 1.5 ~ 2
            )
        }
    )

## Allele counted in UKBB
snp_ea <- data.frame(
    col_id = 2:ncol(snp_count),
    col_n = names(snp_count)[-1]
) %>%
    separate(col_n, c("rsid", "ea_ukbb"), sep = "_", remove = F)

## Checking alignment
aligned_ea <- inner_join(
    snp_ea,
    select(mix, rsid, ea, nea),
    by = "rsid"
) %>%
    filter(ea_ukbb == ea | ea_ukbb == nea) %>%
    mutate(to_change = ea_ukbb != ea)

## Aligning counts and keeping columns with relevant alleles
snp_count <- snp_count %>%
    mutate_at(
        vars(aligned_ea$col_id[aligned_ea$to_change]),
        function(x) 2 - x
    ) %>%
    select(c(1, aligned_ea$col_id))

## Finally, parsing names
names(snp_count) <- gsub("_.*$", "", names(snp_count))

## Saving
save(snp_count, file = "../files/snp_count.RData")
