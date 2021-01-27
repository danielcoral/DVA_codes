#### Calculating scores

## Required libraries
library(tidyverse)

## SNP count - Already aligned to the allele increasing BMI
snp_count <- rio::import("../files/snp_count.RData")

## Concordant and discordant SNPs
mix <- rio::import("../files/mix.txt")

## Summing up allele counts and weights per profile
scores <- lapply(
    ## Codes for profiles
    ## 0 = Concordant, 1 = Discordant
    0:1,
    function(profile){
        dat <- mix[mix$disc == profile,]
        snps <- dat$rsid
        weights <- dat$beta.bmi
        counts <- select(snp_count, all_of(snps))
        res <- apply(counts, 1, sum)
        as.vector(scale(res))
    }
) %>%
    bind_cols() %>%
    `names<-`(c("prs_conc", "prs_disc"))

## Joining and adding identifier
all_scores <- bind_cols(eid = snp_count$IID, scores)

save(all_scores, file = "../files/all_scores.RData")
