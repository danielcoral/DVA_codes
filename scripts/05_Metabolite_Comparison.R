library(readr)
library(dplyr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(purrr)
library(ggplot2)
metab_scan <- read_tsv("../data/metab_scan.tsv", show_col_types = FALSE)
print(head(data.frame(metab_scan)))
profile_meta_metab <- metab_scan %>%
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
head(profile_meta_metab)
write_tsv(profile_meta_metab, "../data/profile_meta_metab.tsv")
fdr_diff_metab <- profile_meta_metab %>%
    mutate(across(c(p_meta, p_comp), ~ p.adjust(.x, "fdr"), .names = "{.col}_adj")) %>%
    filter(p_comp_adj < 0.05) %>%
    group_by(grpid) %>%
    filter(any(p_meta_adj < 0.05)) %>%
    ungroup
head(fdr_diff_metab)
map_ids_metab <- metab_scan %>%
    filter(grpid %in% fdr_diff_metab$grpid) %>%
    select(grpid, trait) %>%
    unique
head(map_ids_metab)
b_res <- Boruta::Boruta(x = snptraitmat_metab[,-c(1:2)], y = factor(snptraitmat_metab$disc), maxRuns = 1000)
imp_df <- Boruta::attStats(b_res)
imp_df$grpid <- rownames(imp_df)
rownames(imp_df) <- NULL
imp_df <- imp_df %>%
    inner_join(map_ids_metab, by = "grpid") %>%
    filter(decision == "Confirmed")
head(imp_df)
imp_df %>%
    write_tsv("../data/imp_df_metab.tsv")
metab_scan %>%
    filter(grpid %in% imp_df$grpid) %>%
    select(id) %>%
    unique %>%
    write_tsv("../data/diff_metabolites.tsv")
