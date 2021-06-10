## Reading VCF files produced

library(tidyverse)
library(vcfR)

## Concordant and discordant SNPs
mix <- rio::import("~/dva/files/mix.tsv")

## SNP clusters
snpclus <- rio::import("~/dva/files/snpclus.tsv")

## VCF files containing genotypes
vcf_files <- list.files("~/dva/files/ukbgeno", full.names = TRUE)

## Extracting allele counts
allele_counts <- map(
    vcf_files,
    ~{
        message(paste("Analysing", .x))
        v <- read.vcfR(.x)
        g <- extract.gt(v, element = 'G', as.numeric = TRUE) %>%
            t %>% data.frame %>% rownames_to_column("eid")
        message("Done!\n")
        return(g)
    }
) %>%
    reduce(inner_join)

rio::export(allele_counts, "~/dva/files/allele_counts.tsv")

allele_counts <- rio::import("~/dva/files/allele_counts.tsv")

## Calculating concordant and discordant scores
scores <- mix %>%
    group_by(disc) %>%
    group_map(~transmute(allele_counts, eid,
                         score = rowSums(across(all_of(.x$rsid)))) %>%
                  `names<-`(c("eid", ifelse(.y$disc == 1, "disc", "conc")))) %>%
    reduce(inner_join)

## Calculating scores for each SNP cluster
snpclus_scores <- snpclus %>%
    mutate(label = gsub(".*\\(|\\)", "", label)) %>%
    group_by(disc, clus_nm) %>%
    mutate(clusid = paste0("snpclus_",
                           ifelse(disc == 1, "disc", "conc"),
                           "_n", n())) %>%
    group_map(~ transmute(allele_counts, eid,
                          score = rowSums(across(all_of(.x$label)))) %>%
                  `names<-`(c("eid", unique(.x$clusid)))) %>%
    reduce(inner_join)

all_scores <- inner_join(scores, snpclus_scores)

rio::export(all_scores, "~/dva/files/all_scores.tsv")
