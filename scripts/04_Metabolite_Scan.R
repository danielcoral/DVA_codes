## Scans for effects of concordant and discordant SNPs on metabolites

library(readr)
library(dplyr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(purrr)
library(ggplot2)
library(ieugwasr)

## Importing concordant and discordant SNPs from previous analysis
mix <- rio::import("../data/mix.tsv")

## And their proxies up to an r2 = 0.5 over 500 kb
## Note: this table include the reference SNPs
proxies <- rio::import("../data/proxies.tsv")

## Batches to use:
batchestouse <- c("met-a", "met-c", "met-d")

## Note 1: All these studies were performed in European populations
## Note 2: All these studies are in continuous traits

## Information of the studies
metab_info <- gwasinfo() %>%
    filter(grepl(paste(batchestouse, collapse = "|"), id))

metab_scan_ref <- lapply(
    split(mix$rsid, ceiling(1:nrow(mix)/10)),
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
    }
) %>%
    bind_rows()

rio::export(metab_scan_ref, "../data/metab_scan_ref.tsv")

metab_scan_ref <- rio::import("../data//metab_scan_ref.tsv")

## Completing missing associations with proxies
## Note: This step is only for proteins and metabolites
metab_scan_prox <- metab_scan_ref %>%
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

rio::export(metab_scan_prox, "../data/metab_scan_prox.tsv")

metab_scan_prox <- rio::import("../data/metab_scan_prox.tsv")

metab_scan <- list(metab_scan_ref, metab_scan_prox) %>%
    bind_rows %>%
    select(-c(chr, n_ld, ref_rsid)) %>%
    inner_join(metab_info) %>%
    inner_join(proxies) %>%
    ## Harmonization - alignment to BMI increasing allele
    mutate(across(c(se, p, position, n, beta, eaf), as.numeric),
           harmon = case_when(ea == ea.bmi & nea == nea.bmi ~ 1,
                              ea == nea.bmi & nea == ea.bmi ~ -1),
           beta = beta * harmon,
           eaf = ifelse(harmon == 1, eaf, 1 - eaf),
           mafstudy = coalesce(ifelse(eaf <= 0.5, eaf, 1 - eaf), maf),
           n = coalesce(n, sample_size)) %>%
    filter(!is.na(beta), abs(maf - mafstudy) < 0.2, n > 500) %>%
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
           se = 1/sqrt(2 * maf * (1 - maf) * (n + (zscore^2))),
           beta = zscore * se)

rio::export(metab_scan, "../data/metab_scan.tsv")

metab_scan <- rio::import("../data/metab_scan.tsv")

