## Difference in protein levels between concordant and discordant profiles

library(tidyverse)
library(ieugwasr)
library(ggrepel)


## Importing concordant and discordant SNPs from previous analysis
mix <- rio::import("../data/mix.tsv")

## And their proxies up to an r2 = 0.5 over 500 kb
proxies <- rio::import("../data/proxies.tsv")

## Batch to use
batchtouse <- "prot-a"

## Information of the studies
prot_info <- gwasinfo() %>%
    filter(grepl(batchtouse, id))

## Associations of reference SNPs
prot_scan_ref <- lapply(
    split(mix$rsid, ceiling(1:nrow(mix)/10)),
    function(snps){
        message("Querying 10 SNPs...", appendLF = FALSE)
        snpres <- phewas(snps, pval = 1, batch = batchtouse) %>%
            mutate(across(.fns = as.character))
        Sys.sleep(1)
        message("\n10 SNPs completed.\n")
        return(snpres)
    }
) %>%
    bind_rows()

rio::export(prot_scan_ref, "../data/prot_scan_ref.tsv")

## Completing missing associations with proxies
prot_scan_prox <- prot_scan_ref %>%
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

rio::export(prot_scan_prox, "../data/prot_scan_prox.tsv")

## Harmonization
prot_scan <- list(prot_scan_ref, prot_scan_prox) %>%
    bind_rows %>%
    select(-c(chr, n_ld, ref_rsid)) %>%
    inner_join(prot_info) %>%
    inner_join(proxies) %>%
    ## Alignment to BMI increasing allele
    mutate(across(c(se, p, position, n, beta, eaf), as.numeric),           
           harmon = case_when(ea == ea.bmi & nea == nea.bmi ~ 1,
                              ea == nea.bmi & nea == ea.bmi ~ -1),
           beta = beta * harmon,
           eaf = ifelse(harmon == 1, eaf, 1 - eaf),
           mafstudy = coalesce(ifelse(eaf <= 0.5, eaf, 1 - eaf), maf),
           n = coalesce(n, sample_size)) %>%
    filter(!is.na(beta), abs(maf - mafstudy) < 0.2) %>%
    group_by(id, ref_rsid) %>%
    arrange(proxy, desc(r2), desc(n), p) %>%
    slice(1) %>%
    ungroup() %>%
    ## Calculating Z-scores and standardizing coefficients
    mutate(zscore = beta / se,
           se = 1/sqrt(2 * maf * (1 - maf) * (n + (zscore^2))),
           beta = zscore * se)

rio::export(prot_scan, "../data/prot_scan.tsv")

prot_scan <- rio::import("../data/prot_scan.tsv")
