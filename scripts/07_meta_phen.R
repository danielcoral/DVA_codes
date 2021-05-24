## Pooling concordant and discordant effects

library(tidyverse)

pheno_scan <- rio::import("~/dva/files/pheno_scan.tsv")

meta_res <- pheno_scan %>%
    group_by(id, disc) %>%
    group_modify(~tryCatch({
        res <- meta::metagen(TE = beta, seTE = se, data = .x,
                             comb.fixed = F, comb.random = T,
                             method.tau = "PM", hakn = F,
                             prediction = F)
        message(paste(paste(.y, collapse = " "), "- done."))
        return(with(res,
                    data.frame(b_meta = TE.random, se_meta = seTE.random,
                               conf_low_meta = lower.random, conf_high_meta = upper.random,
                               p_meta = pval.random, Qv = Q, Qdf = df.Q, Qp = pval.Q,
                               Tau = tau)))
    }, error = function(e) {
        message(paste(paste(.y, collapse = " "), 
                      " could not be completed. Reason:\n", e))
        return(data.frame(NULL))
    })) %>%
    group_by(id) %>%
    filter(dplyr::n() == 2) %>%
    mutate(b_comp = abs(b_meta[1] - b_meta[2]),
           b_comp_se = sqrt(se_meta[1]^2 + se_meta[2]^2),
           p_comp = pnorm(-abs(b_comp/b_comp_se))) %>%
    ungroup()

rio::export(meta_res, "~/dva/files/meta_res.tsv")
