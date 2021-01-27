#### Phenome wide scan

### The tidyverse
library(tidyverse)
library(data.table)

### BMI-T2D intersection
mix <- vroom::vroom("../files/mix.txt")


#### Phenoscanner query ####

phenome_scan <- phenoscanner::phenoscanner(snpquery = mix$rsid, pvalue = 1,
                                           proxies = "EUR", r2 = 0.5)

## Phenoscanner association
phen_res <- phenome_scan$results

fwrite(phen_res, "../files/phen_res.txt")

## Phenoscanner SNP information - including proxies
phen_snps <- phenome_scan$snps

fwrite(phen_snps, "../files/phen_snps.txt")

rm(phenome_scan)

##*********************************
## Harmonizing association results 
##*********************************

## QC of SNPs
qc_snps <- phen_snps %>%
    transmute(
        rsid, a1, a2,
        maf = as.numeric(eur),
        maf = ifelse(maf < 0.5, maf, 1 - maf)
    ) %>%
    filter(!is.na(maf)) %>%
    filter(maf > 0.01) %>%
    filter(stringr::str_count(a1) == 1 & stringr::str_count(a2) == 1) %>%
    filter(
        !(((a1 %in% c("A","T") & a2 %in% c("A","T")) & maf > 0.3) |
          ((a1 %in% c("C","G") & a2 %in% c("C","G")) & maf > 0.3))
    ) %>%
    select(rsid, maf)

fwrite(qc_snps, file = "../files/qc_snps.txt")

premat_phen <- phen_res %>%
    ## 1. Setting column classes and harmonizing trait names
    mutate_at(
        vars(beta, se, p, n, n_cases, n_controls, proxy, r2),
        as.numeric
    ) %>%
    mutate(
        n_min = pmin(n_cases, n_controls),
        trait = tolower(trait)
    ) %>%
    ## 2. Filtering for European ancestry
    filter(ancestry=="European") %>%
    ## 3. Adding SNP information
    inner_join(qc_snps, by = "rsid") %>%
    ## 4. Retaining studies with n > 500 or with at 25 minor alleles
    filter((n_min == 0 & n > 500) | (n_min != 0 & n_min * maf > 25)) %>%
    ## 5. Harmonizing p values and direction
    mutate(
        p = ifelse(p < 1e-300, 1e-300, p),
        direction = case_when(
            direction == "+" ~ 1,
            direction == "-" ~ - 1,
            direction == "0" ~ 0,
            T ~ 2
        )
    ) %>%
    ## 6. Removing columns without direction
    filter(!direction == 2) %>%
    ## 7. Selecting SNP-trait pairs
    arrange(snp, trait, proxy, desc(r2), desc(n)) %>%
    distinct(snp, trait, .keep_all = T) %>%
    ## 8. Only including complete traits with at least 1 FDR 5% hit
    group_by(trait) %>%
    filter(n() == nrow(mix)) %>%
    mutate(p_adj = p.adjust(p, "fdr")) %>%
    filter(sum(p_adj < 0.05) >= 1) %>%
    ungroup() %>%
    ## 9. Calculating z-score
    mutate(zscore = direction * abs(qnorm(p / 2))) %>%
    ## 10. Aligning to the BMI increasing allele
    inner_join(
        mix %>%
        select(rsid, ea, nea),
        by = c("snp" = "rsid")
    ) %>%
    mutate(
        harmon = case_when(
            ref_a1 == ea & ref_a2 == nea ~ 1,
            ref_a1 == nea & ref_a2 == ea ~ 2,
            T ~ 0
        ),
        aligned.z = ifelse(harmon == 1, zscore, -zscore)
    ) %>%
    filter(!harmon == 0) %>%
    ## Setting a code for each trait
    group_by(trait) %>%
    mutate(id = paste0("v", cur_group_id())) %>%
    ungroup()

#### Manually identified traits equivalent to T2D
traits_to_exclude <- c(
    "2 hour fasting glucose",
    "fasting glucose",
    "diabetes diagnosed by doctor",
    "eye problems or disorders: diabetes related eye disease",
    "hba1c",
    "medication for cholesterol, blood pressure or diabetes: insulin",
    "medication for cholesterol, blood pressure or diabetes: none of the above",
    "no treatment with medication for cholesterol, blood pressure, diabetes, or take exogenous hormones",
    "insulin-dependent diabetes mellitus",
    "non-insulin-dependent diabetes mellitus",
    "self-reported diabetes",
    "self-reported diabetic eye disease",
    "self-reported type 2 diabetes",
    "started insulin within one year diagnosis of diabetes",
    "treatment with gliclazide",
    "treatment with insulin",
    "treatment with insulin product",
    "treatment with metformin",
    "treatment with pioglitazone",
    "treatment with rosiglitazone",
    "type ii diabetes",
    "type ii diabetes adjusted for bmi",
    "illnesses of father: diabetes",
    "illnesses of father: none of the above, group 1",
    "illnesses of mother: diabetes",
    "illnesses of mother: none of the above, group 1",
    "illnesses of siblings: diabetes",
    "illnesses of siblings: none of the above, group 1"
)

premat_phen <- premat_phen %>%
    filter(!trait %in% traits_to_exclude)

fwrite(premat_phen, "../files/premat_phen.txt")

#### Setting the identifier ####

id_phen <- premat_phen %>%
    select(id, trait, efo, n_cases) %>%
    arrange(id, desc(efo)) %>%
    distinct(id, .keep_all = T) %>%
    mutate(efo = na_if(efo,"-")) %>%
    arrange(as.numeric(gsub("v","",id)))

fwrite(id_phen, "../files/id_phen.txt")

#### Matrix construction ####
mat_phen <- premat_phen %>%
    select(snp, id, aligned.z) %>%
    pivot_wider(names_from = id, values_from = aligned.z) %>%
    mutate(disc = mix$disc[match(snp, mix$rsid)]) %>%
    select(snp, disc, everything()) %>%
    arrange(desc(disc))

fwrite(mat_phen, file = "../files/mat_phen.txt")

#***********************
# SES distrubution
# Binary vs. Continuous 
#***********************

#### Plot of effect size distribution ####

es_dis <- premat_phen %>%
    transmute(
        `Z-score` = aligned.z,
        type = ifelse(n_cases == 0, "Continuous", "Binary")
    ) %>%
    ggplot(aes(`Z-score`)) +
    geom_density() +
    facet_wrap(~type) +
    ylab("Density")

ggsave("../plots/es_dis.png", es_dis)
