## Blood eQTL effects of concordant and discordant SNPs

library(tidyverse)
library(ieugwasr)
library(ggrepel)

## Importing concordant and discordant SNPs from previous analysis
mix <- rio::import("../data/mix.tsv")

## And their proxies up to an r2 = 0.5 over 500 kb
## Note: this table include reference SNPs
proxies <- rio::import("../data/proxies.tsv")

## For blood eqtl the only batch to use:
batchtouse <- "eqtl-a"

## Information of studies
bloodeqtl_info <- gwasinfo() %>%
    filter(grepl(batchtouse, id))

## Scanning
bloodeqtl_res <- lapply(
    split(proxies$rsid, ceiling(1:nrow(proxies)/10)),
    function(snps){
        message("Querying 10 SNPs...", appendLF = FALSE)
        ## Here we use the significance threshold of the original paper (Vosa et al 2018)
        res <- phewas(snps, pval = 1.83e-5, batch = batchtouse) %>%
            mutate(across(.fns = as.character))
        Sys.sleep(1)
        message("done!")
        return(res)
    }
) %>%
    bind_rows()

rio::export(bloodeqtl_res, "../data/bloodeqtl_res.tsv")

## Harmonizing
bloodeqtl_scan <- bloodeqtl_res %>%
    select(-chr) %>%
    inner_join(bloodeqtl_info) %>%
    inner_join(proxies) %>%
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

## Adding gene names
sig_genes <- unique(bloodeqtl_scan$trait)

ensembl_genes <- biomaRt::useEnsembl(biomart = "ensembl",
                                     dataset = "hsapiens_gene_ensembl")

gene_names <- biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                             filters = "ensembl_gene_id",
                             values = sig_genes,
                             mart = ensembl_genes) %>%
    `names<-`(c("trait", "gene_name"))

## Adding names
bloodeqtl_scan <- bloodeqtl_scan %>%
    left_join(gene_names)

rio::export(bloodeqtl_scan, "../data/bloodeqtl_scan.tsv")

bloodeqtl_scan <- rio::import("~/dva/files/bloodeqtl_scan.tsv")
