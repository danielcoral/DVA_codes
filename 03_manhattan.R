## Manhattan plot

library(tidyverse)
library(hudson)

mix <- rio::import("../files/mix.tsv")

bmi_t2d <- rio::import("../files/bmi_t2d.txt") %>%
    mutate(across(c(p.bmi, p.t2d), function(x) ifelse(x < 1e-300, 1e-300, x))) %>%
    filter(p.bmi < 0.05 | p.t2d < 0.05)

bmi <- select(bmi_t2d, SNP = rsid, CHR = chr, POS = hg19_pos, pvalue = p.bmi)
t2d <- select(bmi_t2d, SNP = rsid, CHR = chr, POS = hg19_pos, pvalue = p.t2d)

gmirror(top = bmi, bottom = t2d,
        tline = 5e-8, bline = 5e-8,
        toptitle = "BMI", bottomtitle = "T2D",
        highlight_snp = mix$rsid, highlighter="green",
        file = "../plots/gmirror.png")


