#### CoVVSURF results - Phenome matrix

library("tidyverse")

premat_phen <- rio::import("../files/premat_phen.txt")   ### The effects 
res_phen <- rio::import("../files/res_phen.RData") ### The results
mix <- rio::import("../files/mix.tsv")

## Clusters of best model
bestmodel <- res_phen$vsel

## Extracting cluster membership 
resvars_phen <- res_phen$ptree$cluster %>%
    data.frame %>%
    rownames_to_column() %>%
    `names<-`(c("id", "cluster")) %>%
    group_by(cluster) %>%
    group_modify(~inner_join(.x,
                             res_phen$ptree$var[.y$cluster][[1]] %>%
                             data.frame %>%
                             rownames_to_column("id"),
                             by = "id")) %>%
    ungroup %>%
    ## Clusters included in best model found
    filter(cluster %in% bestmodel) %>%
    inner_join(premat_phen, by = "id") %>%
    group_by(cluster, squared.loading, trait, disc) %>%
    group_modify(~{
        res <- meta::metagen(TE = beta, seTE = se, data = .x,
                             comb.fixed = F, comb.random = T,
                             method.tau = "PM", hakn = F,
                             prediction = F)
        with(res,
             data.frame(b = TE.random, se = seTE.random,
                        ci_lo = lower.random, ci_up = upper.random,
                        p = pval.random, Qv = Q, Qdf = df.Q, Qp = pval.Q,
                        Tau = tau))
    }) %>%
    group_by(trait) %>%
    mutate(bcomp = abs(b[1] - b[2]), bcomp_se = sqrt(se[1]^2 + se[2]^2),
           pcomp = pnorm(-abs(bcomp/bcomp_se))) %>%
    ungroup %>%
    mutate(p_adj = p.adjust(pcomp, "fdr")) %>%
    arrange(match(cluster, bestmodel), desc(squared.loading))

rio::export(resvars_phen, "../files/resvars_phen.tsv")

## Table to include in supplementary data
resvars_phen %>%
    pivot_wider(c(cluster:trait),
                names_from = disc, values_from = c(b:p_adj)) %>%
    select(c(cluster:trait), ends_with("0"), ends_with("1"),
           -c(contains("comp_0"))) %>%
    mutate(cluster = factor(recode(cluster,
                                   "1" = "Hypertension",
                                   "2" = "Dyslipidaemia",
                                   "4" = "Peripheral fat mass",
                                   "5" = "WHR"),
                            levels = c("Hypertension", "Dyslipidaemia",
                                       "Peripheral fat mass", "WHR"))) %>%
    arrange(cluster, desc(squared.loading)) %>%
    rio::export("../files/clusters.xlsx")
