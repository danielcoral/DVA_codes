## Building SNP-trait matrix using results from subtype specific meta-analysis

library(tidyverse)

pheno_scan <- rio::import("~/dva/files/pheno_scan.tsv")

meta_res <- rio::import("~/dva/files/meta_res.tsv")

premat_phen <- pheno_scan %>%
    inner_join(meta_res) %>%
    {
        ids_to_keep <- group_by(., tolower(trait), category, sex) %>%
            arrange(desc(sample_size)) %>%
            slice(1) %>%
            pull(id)
        filter(., id %in% ids_to_keep)
    } %>%        
    mutate(across(c(p_meta, p_comp), ~ p.adjust(.x, "fdr"),
                  .names = "{.col}_adj")) %>%
    filter(p_comp_adj < 0.1) %>%
    group_by(id) %>%
    filter(any(p_meta_adj < 0.1)) %>%
    ungroup

rio::export(premat_phen, "~/dva/files/premat_phen.tsv")

mat_phen <- premat_phen %>%
    transmute(ref_rsid, disc, grpid, zscore) %>%
    pivot_wider(names_from = grpid, values_from = zscore) %>%
    arrange(desc(disc))

dim(mat_phen)

rio::export(mat_phen, "~/dva/files/mat_phen.tsv")
