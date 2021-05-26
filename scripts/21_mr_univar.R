## Performing MR

library(tidyverse)
library(MendelianRandomization)

mrdat_all <- rio::import("~/dva/files/mrdat_all.tsv")
mrdat_bmisig <- rio::import("~/dva/files/mrdat_bmisig.tsv")

## Custom MR Egger function
## (Since preset functions change the orientation of variables)
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

mr_analysis <- function(dat,...){
    ## Main: IVW - Evidence of an effect of the exposure on the outcome:
    ivw_res <- dat %>%
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
    egger_res <- dat %>%
        group_by(exposure) %>%
        group_modify(~ eggerfx(.x$beta.exposure, .x$beta.outcome,
                               .x$se.exposure, .x$se.outcome))
    ## Weighted median: Majority of instruments are valid
    wmedian_res <- dat %>%
        TwoSampleMR::mr(method_list = "mr_weighted_median") %>%
        mutate(method = "wmedian",
               ci_lo = b - qnorm(1-0.05/2)*se,
               ci_up = b + qnorm(1-0.05/2)*se) %>%
        select(-c(id.exposure, id.outcome, outcome))
    ## ConMix: Plurality of variants are valid
    conmix_res <- dat %>%
        group_by(exposure) %>%
        group_modify(~{
            lapply(seq(1,2,.1),
                   function(psi_val){
                       tryCatch({
                           res <- mr_input(bx = .x$beta.exposure, bxse = .x$se.exposure,
                                           by = .x$beta.outcome, byse = .x$se.outcome) %>%
                               mr_conmix(psi = psi_val * sd(.x$beta.outcome/.x$beta.exposure))
                           data.frame(b = res@Estimate,
                                      ci_lo = res@CILower, ci_up = res@CIUpper,
                                      pval = res@Pvalue, nsnp = res@SNPs, psi = psi_val,
                                      psitrue = res@Psi, method = "conmix") %>%
                               mutate(anymulti = na_if(any(duplicated(psi)), FALSE))
                       }, error = function(e){
                           message(paste(e, "\nHappened during ConMix analysis of",
                                         .y$exposure, "- PSI of", psi_val))
                           return(NULL)
                       })
                   }) %>%
                bind_rows()
        }) %>%
        ungroup()
    ## Mean F-statistic
    Fstat <- dat %>%
        group_by(exposure) %>%
        summarise(f_mean = mean((beta.exposure^2) / (se.exposure^2)),
                  .groups = "drop")
    ## Isq_GX
    isq_gx <- dat %>%
        group_by(exposure) %>%
        summarise(isq = TwoSampleMR::Isq(beta.exposure, se.exposure),
                  .groups = "drop")
    ## Joining results
    mr_res <- bind_rows(ivw_res, egger_res, wmedian_res) %>%
        inner_join(Fstat) %>% inner_join(isq_gx)
    return(mr_res)
}

mr_res_univariate <- map_dfr(list(all = mrdat_all,
                                  bmi_sig = mrdat_bmisig),
                             mr_analysis, .id = "analysis")

rio::export(mr_res_univariate, "~/dva/files/mr_res_univariate.tsv")

mr_profiles <- mrdat_bmisig %>%
    mutate(disc = ifelse(beta.outcome < 0, "Discordant", "Concordant")) %>%
    group_by(disc) %>%
    group_modify(mr_analysis)

rio::export(mr_profiles, "~/dva/files/mr_profiles.tsv")
