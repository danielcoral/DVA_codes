## Phenome-wide scan of concordant and discordant SNPs

library(tidyverse)
library(ieugwasr)
library(genetics.binaRies)
library(gwasrapidd)
library(easyPubMed)

##********************************##
## Concordant and discordant SNPs ##
##********************************##
mix <- rio::import("~/dva/files/mix.tsv")

##*********##
## Proxies ##
##*********##

## Files to write PLINK output
out <- "~/dva/files/mix_prox"
targetsname <- paste0(out, ".targets")
outname <- paste0(out, ".targets.ld.gz")

## SNPs to query
utils::write.table(mix$rsid, file=targetsname,
                   row=FALSE, col=FALSE, qu=FALSE)

## PLINK command
cmd <- paste0(
    genetics.binaRies::get_plink_binary(),
    " --bfile ~/dva/files/1kg_ref/EUR", ## Reference
    " --r2 in-phase with-freqs gz",  ## R2 with alleles(PHASE) and MAF
    " --ld-snp-list ", targetsname,  ## SNPs to query
    " --ld-window-kb 500",          ## 500 KB window
    " --ld-window 999999",           ## N of SNPs between each pair
    " --out ", targetsname
)

system(cmd)

## Generating proxies table
proxies <- readr::read_table2(pipe(paste("gunzip -c ", outname))) %>%
    select(c(3, 7, 8, 9, 10)) %>%
    `names<-`(c("ref_rsid", "rsid", "phase", "maf", "r2")) %>%
    filter(r2 > 0.5) %>%                     ## r2 threshold
    filter(maf > 0.01) %>%                   ## MAF threshold
    mutate(phase = gsub("/", "", phase)) %>% ## Remove INDELs
    filter(nchar(phase) == 4) %>%
    ## Formatting alleles (a is reference, b is proxy)
    separate(phase, into = paste0(rep(c("a", "b"), 2),
                                  rep(c("1", "2"), each = 2)),
             sep = c(1,2,3)) %>%
    ## Identifying reference and proxies
    mutate(proxy = ref_rsid != rsid,
           ## Palindromic SNPs
           pal1 = b1 %in% c("A", "T") & b2 %in% c("A", "T"),
           pal2 = b1 %in% c("C", "G") & b2 %in% c("C", "G"),
           ambig = (pal1 | pal2) & maf > 0.3) %>%
    filter(!ambig) %>%
    ## Aligning to the BMI increasing allele
    inner_join(mix, by = c("ref_rsid" = "rsid")) %>%
    mutate(ea.bmi = case_when(a1 == ea & a2 == nea ~ b1,
                              a1 == nea & a2 == ea ~ b2),
           nea.bmi = case_when(a1 == ea & a2 == nea ~ b2,
                               a1 == nea & a2 == ea ~ b1)) %>%
    filter(!is.na(ea.bmi), !is.na(nea.bmi)) %>%
    select(-c(a1,b1,a2,b2,ea,nea,ambig,pal1,pal2)) %>%
    ## Number of proxies per lead SNP
    group_by(ref_rsid) %>%
    mutate(n_ld = dplyr::n()) %>%
    ungroup %>%
    arrange(n_ld)

rio::export(proxies, "~/dva/files/proxies.tsv")

proxies <- rio::import("~/dva/files/proxies.tsv")

##********##
## Traits ##
##********##

## MRC IEU database is divided in batches:
batchids <- batches()

## Examining batches
batchids %>%
    transmute(id,
              description = str_trunc(description, 60))

## Selecting batches
batchestouse <- c("ebi-a", ## Exported from EBI
                  "ieu-a", "ieu-b", ## Harmonized by IEU
                  "ukb-a", "ukb-b", "ukb-d") ## UK Biobank analysis

## Obtaining metadata of traits
availgwas <- gwasinfo()

## Metadata of biomarkers (from Neale lab)
## From this we extract sample size (absent in trait metadata)
biomarker_url <- paste0("https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/",
                        "round2/annotations/biomarkers.both_sexes.tsv.bgz")
biomarkerdat <- read_delim(pipe(paste0("curl -s ", biomarker_url, " | gunzip -c")),
                           delim = "\t") %>%
    transmute(id = paste0("ukb-d-", phenotype), biomarker_n = n_non_missing)

traitinfo <- availgwas %>%
    ## Only for the batches selected
    filter(str_detect(id, paste(batchestouse, collapse = "|")),
           ## Studies performed in European populations
           population == "European" | 
               ## Of population == NA, only Wojcik et al. is not a study in Europeans
               (population == "NA" & pmid != 31217584) |
               ## Mixed population may include European descent
               ## Only including those for which we can check this
               (population == "Mixed" & !is.na(pmid)),
           ## Also removing categorical variables
           ## (Reason: No data of their sample sizes within categories)
           category != "Categorical Ordered") %>%
    left_join(biomarkerdat) %>%
    mutate(sample_size = coalesce(sample_size,
                                  ncase + ncontrol,
                                  biomarker_n),
           category = case_when(category %in% c("Binary", "Disease") ~ "Binary",
                                category %in% c("Risk factor", "NA") & 
                                    !is.na(ncase) & !is.na(ncontrol) ~ "Binary",
                                TRUE ~ "Continuous"),
           sex = case_when(sex == "Males" ~ "M", sex == "Females" ~ "F",
                           TRUE ~ "Both"),
           n_min = pmin(ncase, ncontrol)) %>%
    filter((category == "Continuous" & sample_size > 1000) |
               (category == "Binary" & sample_size > 1000 & n_min > 100)) %>%
    select(id, trait, year, consortium, author, sex, pmid, population, unit,
           nsnp, sample_size, category, ncase, ncontrol, n_min)

## Registering percentage of European descent in studies with population == Mixed

## PubMed IDs of studies with mixed populations
popmixed_pmid <- traitinfo %>%
    filter(population == "Mixed") %>%
    pull(pmid) %>%
    unique

## Extracting metadata for these studies from GWAS Catalog
popmixed_info <- get_studies(pubmed_id = as.character(unique(popmixed_pmid)))

## Shaping metadata
popmixed_tab <- popmixed_info@publications %>%
    select(pubmed_id, study_id) %>%
    inner_join(popmixed_info@ancestral_groups) %>%
    inner_join(popmixed_info@ancestries) %>%
    inner_join(select(popmixed_info@studies, study_id, reported_trait)) %>%
    transmute(pmid = pubmed_id, study_id, reported_trait,
              ancestry = replace_na(ancestral_group, "unspecified"),
              ancestry_n = number_of_individuals) %>%
    group_by(pmid, study_id, reported_trait) %>%
    arrange(ancestry != "European") %>%
    distinct(ancestry_n, .keep_all = TRUE) %>%
    group_by(ancestry, .add = TRUE) %>%
    summarise(ancestry_total_n = sum(ancestry_n)) %>%
    summarise(n_eur_study = ancestry_total_n[ancestry == "European"],
              n_total = sum(ancestry_total_n),
              istotal = dplyr::n() == 1) %>%
    mutate(perc_eur = n_eur_study / n_total) %>%
    arrange(perc_eur)

## Those not found in GWAS catalog can be manually checked
## We do this using `easyPubMed`
pmid_add <- setdiff(popmixed_pmid, popmixed_tab$pmid)

pmid_add

lapply(pmid_add, 
       function(x){
           res <- get_pubmed_ids(x)
           fetch_pubmed_data(res, format = "abstract", retmax = 1)
       })

## Then we can manually extract the proportion of sample that is from European descent
## In all, the majority of the sample analysed is from European descent
popmixed_add <- data.frame(pmid = pmid_add,
                           n_eur_study = c(13435+29269,
                                           27206+57574,
                                           63746+130681-(6557+3625+4789),
                                           2668),
                           n_total = c(13435+29269+2385+2193,
                                       27206+57574,
                                       63746+130681+6557+3625+4789,
                                       3829),
                           istotal = c(FALSE, TRUE, FALSE, FALSE)) %>%
    mutate(perc_eur = n_eur_study / n_total) %>%
    inner_join(select(traitinfo, pmid, reported_trait = trait))

popmixed_tabf <- bind_rows(popmixed_tab, popmixed_add)

rio::export(popmixed_tabf, "~/dva/files/popmixed_tabf.tsv")

## Filtering traits that contain or are used for diabetes diagnosis
traitinfo %>%
    pull(trait) %>%
    unique %>%
    paste(collapse = " || ")

toremove <- c("Type 2 diabetes", ## T2D diagnosis
              "Type 2 diabetes (adjusted for BMI)",
              "Non-cancer illness code  self-reported: type 2 diabetes",
              "Non-cancer illness code, self-reported: type 2 diabetes",
              "Non-cancer illness code, self-reported: diabetes",
              "Non-cancer illness code  self-reported: diabetes",
              "Diabetes diagnosed by doctor", 
              "Age diabetes diagnosed",
              "Started insulin within one year diagnosis of diabetes",
              "Diagnoses - secondary ICD10: E11.9 Without complications",
              "Diagnoses - secondary ICD10: E10.9 Without complications",
              
              ## Glucose measurements
              "Glucose", "Fasting glucose", "Fasting blood glucose", "Two-hour glucose challenge",
              "Fasting blood glucose adjusted for BMI", 
              "HbA1C", "Glycated haemoglobin", "Glycated hemoglobin levels",
              
              ## Family history of T2D
              "Illnesses of mother: Diabetes",
              "Illnesses of father: Diabetes",
              "Illnesses of siblings: Diabetes",
              "Illnesses of adopted siblings: Diabetes",
              ## Group 1 includes diabetes
              "Illnesses of mother: None of the above (group 1)",
              "Illnesses of father: None of the above (group 1)",
              "Illnesses of siblings: None of the above (group 1)",
              "Illnesses of adopted siblings: None of the above (group 1)",
              
              ## T2D medication
              "Medication for cholesterol  blood pressure or diabetes: Insulin",
              "Medication for cholesterol, blood pressure or diabetes: Insulin",
              "Medication for cholesterol  blood pressure or diabetes: None of the above",
              "Medication for cholesterol, blood pressure or diabetes: None of the above",
              "Medication for cholesterol  blood pressure, diabetes, or take exogenous hormones: Insulin",
              "Medication for cholesterol, blood pressure, diabetes, or take exogenous hormones: Insulin",
              "Medication for cholesterol  blood pressure  diabetes  or take exogenous hormones: None of the above",
              "Medication for cholesterol, blood pressure, diabetes, or take exogenous hormones: None of the above",
              "Treatment/medication code: insulin product",
              "Treatment/medication code: metformin",
              "Treatment/medication code: rosiglitazone",
              "Treatment/medication code: gliclazide",
              
              ## T2D complications (include T2D diagnosis)
              "Eye problems/disorders: Diabetes related eye disease",
              "Glomerular filtration rate in diabetics (creatinine)",
              "Glomerular filtration rate in non diabetics (creatinine)")

traitstouse <- traitinfo %>%
    filter(!trait %in% toremove)

rio::export(traitstouse, "~/dva/files/traitstouse.tsv")

traitstouse <- rio::import("~/dva/files/traitstouse.tsv")

##**********************************##
## Phenome-scan with reference SNPs ##
##**********************************##

pheno_scan_ref <- lapply(split(mix$rsid, ceiling(1:nrow(mix)/10)),
                         function(snps){
                             message("Querying 10 SNPs...", appendLF = FALSE)
                             snpres <- lapply(batchestouse,
                                              function(batchid){
                                                  message("\n+--Querying batch...", appendLF = FALSE)
                                                  batchres <- phewas(snps, pval = 1, batch = batchid) %>%
                                                      mutate(across(.fns = as.character))
                                                  message("batch completed.", appendLF = FALSE)
                                                  Sys.sleep(1)
                                                  return(batchres)
                                              }) %>%
                                 bind_rows()
                             message("\n10 SNPs completed.\n")
                             return(snpres)
                         }) %>%
    bind_rows()

rio::export(pheno_scan_ref, "~/dva/files/pheno_scan_ref.tsv")

pheno_scan_ref <- rio::import("~/dva/files/pheno_scan_ref.tsv")

##**********************************************##
## Completing missing associations with proxies ##
##**********************************************##

pheno_scan_prox <- pheno_scan_ref %>%
    filter(id %in% traitstouse$id) %>%
    group_by(id) %>%
    filter(dplyr::n() < nrow(mix)) %>%
    summarise(ref_rsid = mix$rsid[!mix$rsid %in% rsid],
              .groups = "drop") %>%
    inner_join(proxies) %>%
    filter(proxy) %>%
    group_by(id) %>%
    group_modify(
        function(trait_dat,...){
            message(paste0("\nQuerying phenotype ", unique(trait_dat$id), "..."), 
                    appendLF = FALSE)
            phenores <- tryCatch({
                trait_dat %>%
                    group_by(n_ld, ref_rsid) %>%
                    group_modify(
                        function(prox,...){
                            message(paste0("\n+--Querying reference SNP ",
                                           unique(prox$ref_rsid), "..."), appendLF = FALSE)
                            resprox <- associations(prox$rsid, unique(prox$id), proxies = 0) %>%
                                mutate(across(.fns = as.character)) %>%
                                select(-any_of("id"))
                            if(nrow(resprox) == 0){
                                stop("One of the variants was not found.")
                            } else {
                                message("completed.", appendLF = FALSE)
                                return(resprox)
                            }
                        }, .keep = TRUE
                    )
                }, error = function(e) {
                    message(e, appendLF = FALSE)
                    return(data.frame(NULL))
                })
            message("\n...phenotype completed.")
            return(phenores)
        }, .keep = TRUE
    )

rio::export(pheno_scan_prox, "~/dva/files/pheno_scan_prox.tsv")

pheno_scan_prox <- rio::import("~/dva/files/pheno_scan_prox.tsv")

##*****************##
## Joining results ##
##*****************##

pheno_scan <- list(pheno_scan_ref,
                   pheno_scan_prox) %>%
    bind_rows %>%
    select(-c(chr, n_ld, ref_rsid)) %>%
    ## Adding metadata of proxies and traits
    inner_join(traitstouse) %>%
    inner_join(proxies) %>%
    ## Harmonization - alignment to BMI increasing allele
    mutate(across(c(se, p, position, n, beta), as.numeric),
           harmon = case_when(ea == ea.bmi & nea == nea.bmi ~ 1,
                              ea == nea.bmi & nea == ea.bmi ~ -1),
           beta = beta * harmon,
           eaf = ifelse(harmon == 1, eaf, 1 - eaf),
           mafstudy = coalesce(ifelse(eaf <= 0.5, eaf, 1 - eaf), maf),
           ma_n = 2 * maf * n_min,
           n = coalesce(n, sample_size)) %>%
    filter(!is.na(beta), se > 0, p >= 0, p < 1,
           ## Threshold for valid associations for binary traits
           category == "Continuous" |
               category == "Binary" & ma_n > 25,
           ## Difference in MAF with 1KG EUR < 20%
           abs(maf - mafstudy) < 0.2) %>%
    ## One association per reference SNP - trait pair
    group_by(id, ref_rsid) %>%
    arrange(proxy, desc(r2), desc(n), p) %>%
    slice(1) %>%
    ungroup(ref_rsid) %>%
    ## Only traits with information for all reference SNPs
    filter(length(unique(ref_rsid)) == nrow(mix)) %>%
    mutate(grpid = paste0("v", cur_group_id())) %>%
    ungroup() %>%
    ## Calculating Z-scores and standardizing coefficients
    mutate(zscore = beta / se,
           se = ifelse(category == "Continuous",
                       1/sqrt(2 * maf * (1 - maf) * (n + (zscore^2))), se),
           beta = ifelse(category == "Continuous",
                         zscore * se, beta))
##           beta = ifelse(category == "Continuous",
##                         sqrt((zscore^2/(zscore^2+n-2)) / (2 * maf * (1-maf))) * sign(zscore),
##                         beta),
##           se = beta / zscore)

rio::export(pheno_scan, "~/dva/files/pheno_scan.tsv")
