library(readr)
library(dplyr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(purrr)
library(ggplot2)
prot_scan <- read_tsv("../data/prot_scan.tsv", show_col_types = FALSE)
print(head(data.frame(prot_scan)))
profile_meta_prot <- prot_scan %>%
    ## For a trait, selecting the study with the highest sample size
    mutate(trait = tolower(trait)) %>%
    group_by(trait, category, sex) %>%
    filter(id == id[which.max(sample_size)]) %>%
    ## Meta-analysing SNP effects for each genetic profile in every trait
    group_by(id, disc) %>%
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
    group_by(id) %>%
    filter(dplyr::n() == 2) %>%
    mutate(b_comp = abs(b_meta[1] - b_meta[2]),
           b_comp_se = sqrt(se_meta[1]^2 + se_meta[2]^2),
           p_comp = 2 * pnorm(-abs(b_comp/b_comp_se))) %>%
    ungroup
head(profile_meta_prot)
write_tsv(profile_meta_prot, "../data/profile_meta_prot.tsv")
fdr_diff_prot <- profile_meta_prot %>%
    mutate(across(c(p_meta, p_comp), ~ p.adjust(.x, "fdr"), .names = "{.col}_adj")) %>%
    filter(p_comp_adj < 0.05) %>%
    group_by(id) %>%
    filter(any(p_meta_adj < 0.05)) %>%
    ungroup
head(fdr_diff_prot)
unique(prot_scan$trait[prot_scan$id == unique(fdr_diff_prot$id)])
p1_mr <- prot_scan %>%
    filter(p < 5e-8, disc == 1)
p1_mr %>%
    data.frame %>%
    print
