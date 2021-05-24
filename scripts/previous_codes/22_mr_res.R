## Performing MR

library(tidyverse)
library(MendelianRandomization)

pheno_exp <- rio::import("../files/pheno_exp.tsv")

## Main: IVW
## Evidence of an effect of the exposure on the outcome:
ivw_res <- pheno_exp %>%
    {list(TwoSampleMR::mr(., method_list = "mr_ivw"),
          TwoSampleMR::mr_heterogeneity(., method_list = "mr_ivw"))} %>%
    reduce(inner_join) %>%
    mutate(method = "ivw",
           ci_lo = b - qnorm(1-0.05/2)*se,
           ci_up = b + qnorm(1-0.05/2)*se,
           padj = p.adjust(pval, "fdr")) %>%
    select(-c(id.exposure, id.outcome, outcome))

## MR-Egger: balanced pleiotropy (InSIDE)
## Also tests a dose-response relationship
## Since preset functions change the orientation of variables:
eggerfx <- function(b_exp, b_out, se_exp, se_out){
    mod <- lm(b_out ~ b_exp, weights=1/se_out^2)
    smod <- summary(mod)
    data.frame(method = c("intercept", "egger"),
               b = coefficients(smod)[,1],
               se = coefficients(smod)[,2] / min(1,smod$sigma),
               Q = smod$sigma^2 * (length(b_exp) - 2),
               Q_df = length(b_exp) - 2,
               nsnp = length(b_exp)) %>%
        mutate(pval = 2 * pt(abs(b / se), length(b_exp) - 2, lower.tail = FALSE),
               Q_pval = pchisq(Q, Q_df, lower.tail=FALSE),
               ci_lo = b - qnorm(1-0.05/2)*se,
               ci_up = b + qnorm(1-0.05/2)*se)
}

egger_res <- pheno_exp %>%
    group_by(exposure) %>%
    group_modify(~ eggerfx(.x$beta.exposure, .x$beta.outcome,
                           .x$se.exposure, .x$se.outcome))

## Weighted median: Majority of instruments are valid
wmedian_res <- pheno_exp %>%
    TwoSampleMR::mr(method_list = "mr_weighted_median") %>%
    mutate(method = "wmedian",
           ci_lo = b - qnorm(1-0.05/2)*se,
           ci_up = b + qnorm(1-0.05/2)*se) %>%
    select(-c(id.exposure, id.outcome, outcome))

## ConMix: Plurality of variants are valid
conmix_res <- pheno_exp %>%
    group_by(exposure) %>%
    group_modify(~{
        lapply(seq(1,2,.1),
               function(psi){
                   res <- mr_input(bx = .x$beta.exposure, bxse = .x$se.exposure,
                                   by = .x$beta.outcome, byse = .x$se.outcome) %>%
                       mr_conmix(psi = psi * sd(.x$beta.outcome/.x$beta.exposure))
                   data.frame(b = res@Estimate, ci_lo = res@CILower, ci_up = res@CIUpper,
                              pval = res@Pvalue, nsnp = res@SNPs, psi = psi,
                              psitrue = res@Psi, method = "conmix") %>%
                       mutate(anymulti = na_if(any(duplicated(psi)), FALSE))
               }) %>%
            bind_rows()
    }) %>%
    ungroup()

## Mean F-statistic
Fstat <- pheno_exp %>%
    group_by(exposure) %>%
    summarise(f_mean = mean((beta.exposure^2) / (se.exposure^2)),
              .groups = "drop")

## Isq_GX
isq_gx <- pheno_exp %>%
    group_by(exposure) %>%
    summarise(isq = TwoSampleMR::Isq(beta.exposure, se.exposure),
              .groups = "drop")

## Joining results
mr_res <- bind_rows(ivw_res, egger_res, wmedian_res, conmix_res) %>%
    ## Applying Rucker's framework
    group_by(exposure) %>%
    mutate(pleio.test = any(method == "intercept" & pval < 0.05),
           Qdiff = Q[method == "egger"] - Q[method == "ivw"],
           Qdiff.p = pchisq(Qdiff, df = 1, lower.tail = FALSE)) %>%
    ungroup %>%
    inner_join(Fstat) %>% inner_join(isq_gx)

rio::export(mr_res, "../files/mr_res.tsv")
