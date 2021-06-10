## Multivariable MR of SBP adjusted for fasting insulin

library(tidyverse)
library(ieugwasr)

## MR data
harmondat_all <- rio::import("~/dva/files/harmondat_all.tsv")

## SBP
sbp <- harmondat_all %>%
    filter(trait_short == "SBP")

## Fasting insulin from Manning et al 2012
fi <- associations(sbp$rsid, "ieu-b-115", proxies = 0)

fi <- fi %>%
    select(-c(chr, position, id, trait, eaf)) %>%
    rename_with(~paste(.x, "fi", sep = "."), -rsid)

mvmrdat <- sbp %>%
    select(rsid, ea, nea, MAF,
           beta.sbp = beta.exposure, se.sbp = se.exposure,
           p.sbp = p.exposure,
           beta.bmi, se.bmi, p.bmi,
           beta.t2d = beta.outcome, se.t2d = se.outcome, p.t2d = pval.outcome) %>%
    inner_join(fi) %>%
    mutate(harmon = case_when(ea == ea.fi & nea == nea.fi ~ 1,
                              ea == nea.fi & nea == ea.fi ~ -1,
                              TRUE ~ 0),
           beta.fi = beta.fi * harmon,
           z.fi = beta.fi / se.fi,
           se.fi = 1/sqrt(2 * MAF * (1 - MAF) * (n.fi + (z.fi^2))),
           beta.fi = z.fi * se.fi) %>%
    select(-c(z.fi, harmon, ea.fi, nea.fi, n.fi))

## Clumping
toclump <- select(mvmrdat, rsid, pval = p.sbp)

clumped <- ld_clump(toclump, clump_kb = 500, clump_r2 = 0.01,
                    plink_bin = genetics.binaRies::get_plink_binary(),
                    bfile = "/ludc/Home/daniel_c/dva/files/1kg_ref/EUR")

mvmrdat_c <- filter(mvmrdat, rsid %in% clumped$rsid)

## MVMR functions
source(list.files(pattern = "mrfx.R$"))

mvmr_res <- with(mvmrdat_c, mvmr_analysis(cbind(beta.sbp, beta.fi), beta.t2d,
                                          cbind(se.sbp, se.fi), se.t2d,
                                          exposures = c("sbp", "fi"))) %>%
    tibble

mvmr_res

rio::export(mvmr_res, "~/dva/files/mvmr_res.tsv")
