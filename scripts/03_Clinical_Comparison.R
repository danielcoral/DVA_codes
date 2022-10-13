library(readr)
library(dplyr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(ggplot2)
library(patchwork)
pheno_scan <- read_tsv("../data/pheno_scan.tsv", show_col_types = FALSE)
print(head(data.frame(pheno_scan)))
profile_meta_pheno <- pheno_scan %>%
    ## For a trait, selecting the study with the highest sample size
    mutate(trait = tolower(trait)) %>%
    group_by(trait, category, sex) %>%
    filter(grpid == grpid[which.max(sample_size)]) %>%
    ## Meta-analysing SNP effects for each genetic profile in every trait
    group_by(grpid, disc) %>%
    group_modify(~tryCatch({
        res <- meta::metagen(TE = beta, seTE = se, data = .x,
                             fixed = FALSE, random = TRUE,
                             method.tau = "PM",
                             prediction = FALSE)
        return(with(res,
                    data.frame(b_meta = TE.random, se_meta = seTE.random,
                               conf_low_meta = lower.random, conf_high_meta = upper.random,
                               p_meta = pval.random, Qv = Q, Qdf = df.Q, Qp = pval.Q,
                               Tau = tau)))
    }, error = function(e) {
        return(data.frame(NULL))
    })) %>%
    group_by(grpid) %>%
    filter(dplyr::n() == 2) %>%
    mutate(b_comp = abs(b_meta[1] - b_meta[2]),
           b_comp_se = sqrt(se_meta[1]^2 + se_meta[2]^2),
           p_comp = 2 * pnorm(-abs(b_comp/b_comp_se))) %>%
    ungroup
head(profile_meta_pheno)
write_tsv(profile_meta_pheno, "../data/profile_meta_pheno.tsv")
fdr_diff_pheno <- profile_meta_pheno %>%
    mutate(across(c(p_meta, p_comp), ~ p.adjust(.x, "fdr"), .names = "{.col}_adj")) %>%
    filter(p_comp_adj < 0.05) %>%
    group_by(grpid) %>%
    filter(any(p_meta_adj < 0.05)) %>%
    ungroup
head(fdr_diff_pheno)
write_tsv(fdr_diff_pheno, "../data/fdr_diff_pheno.tsv")
snptraitmat_phen <- pheno_scan %>%
    filter(grpid %in% fdr_diff_pheno$grpid) %>%
    transmute(ref_rsid, disc, grpid, zscore) %>%
    pivot_wider(names_from = grpid, values_from = zscore) %>%
    arrange(desc(disc))
dim(snptraitmat_phen)
map_ids_pheno <- pheno_scan %>%
    filter(grpid %in% fdr_diff_pheno$grpid) %>%
    select(grpid, trait, sex, category) %>%
    unique
head(map_ids_pheno)
b_res <- Boruta::Boruta(x = snptraitmat_phen[,-c(1:2)], y = factor(snptraitmat_phen$disc), maxRuns = 1000)
imp_df <- Boruta::attStats(b_res)
imp_df$grpid <- rownames(imp_df)
rownames(imp_df) <- NULL
imp_df <- imp_df %>%
    inner_join(map_ids_pheno, by = "grpid") %>%
    filter(decision == "Confirmed") %>%
    mutate(trait = ifelse(sex == "Both", trait, paste(trait, sex, sep =" - "))) %>%
    select(-sex)
head(imp_df)
write_tsv(imp_df, "../data/imp_df.tsv")
