##Building PRSs

library(data.table)

## Using 1KG as map

ref <- fread("../files/1kg_ref/EUR.bim", select = c(2,1,4),
             col.names = c("rsid", "chr", "hg19_pos"))

## BMI GWAS from GIANT

bmi_giant <- fread("../files/bmi_giant.txt") 

bminoukb <- bmi_giant[
    p < 5e-8
][
    ref, on = c(SNP = "rsid"), nomatch = 0
][
   ,.(trait = "bmi", rsid = SNP, chr, hg19_pos,
      ea = ifelse(b > 0, A1, A2), pval = p)
]

## BP GWAS from Evangelou et al. (ICBP + UKB)

sbp <- fread("../files/mr_ss/sbp.txt.gz")

dbp <- fread("../files/mr_ss/dbp.txt.gz")

## BP GWAS ICBP only

sbp_icbp <- fread("../files/sbp_icbp.tsv.gz")

dbp_icbp <- fread("../files/sbp_icbp.tsv.gz")

## Only selecting variants significant in ICBP

sbpnoukb <- sbp_icbp[
    p_value < 5e-8
][
    ref, on = c(hm_rsid = "rsid"), nomatch = 0
][
  , MarkerName := paste(chr, hg19_pos, "SNP", sep = ":")
][
    sbp, on = "MarkerName", nomatch = 0
][
  , .(trait = "sbp", rsid = hm_rsid, chr, hg19_pos,
      ea = toupper(ifelse(Effect > 0, Allele1, Allele2)),
      pval = p_value)
]

dbpnoukb <- dbp_icbp[
    p_value < 5e-8
][
    ref, on = c(hm_rsid = "rsid"), nomatch = 0
][
  , MarkerName := paste(chr, hg19_pos, "SNP", sep = ":")
][
    dbp, on = "MarkerName", nomatch = 0
][
  , .(trait = "dbp", rsid = hm_rsid, chr, hg19_pos,
      ea = toupper(ifelse(Effect > 0, Allele1, Allele2)),
      pval = p_value)
]

## Cholesterol from GLGC

hdl_glgc <- fread("../files/hdl_glgc.txt.gz")

ldl_glgc <- fread("../files/ldl_glgc.txt.gz")

tg_glgc <- fread("../files/tg_glgc.txt.gz")

hdlnoukb <- hdl_glgc[, pval := as.numeric(`P-value`)][pval < 5e-8][
    , .(trait = "hdl", rsid,
        chr = as.numeric(gsub("chr|:[0-9]+", "", SNP_hg19)),
        hg19_pos = as.numeric(gsub("chr[0-9]+:", "", SNP_hg19)),
        ea = toupper(ifelse(beta > 0, A1, A2)), pval)
]

ldlnoukb <- ldl_glgc[, pval := as.numeric(`P-value`)][pval < 5e-8][
    , .(trait = "ldl", rsid,
        chr = as.numeric(gsub("chr|:[0-9]+", "", SNP_hg19)),
        hg19_pos = as.numeric(gsub("chr[0-9]+:", "", SNP_hg19)),
        ea = toupper(ifelse(beta > 0, A1, A2)), pval)
]

tgnoukb <- tg_glgc[, pval := as.numeric(`P-value`)][pval < 5e-8][
  , .(trait = "tg", rsid,
      chr = as.numeric(gsub("chr|:[0-9]+", "", SNP_hg19)),
      hg19_pos = as.numeric(gsub("chr[0-9]+:", "", SNP_hg19)),
      ea = toupper(ifelse(beta > 0, A1, A2)), pval)
]

## WHR from GIANT
whr_giant <- fread("../files/whr_giant.txt")

whrnoukb <- whr_giant[
    p < 5e-8
][
   ,.(trait = "whr", rsid = MarkerName, chr = Chr, hg19_pos = Pos,
      ea = ifelse(b > 0, Allele1, Allele2), pval = p)
]

## Liver function from Chambers et al & Prins et al

ggt <- fread("../files/ggt.txt.gz")

ggtnoukb <- ggt[
    P_value < 5e-8
][
    ref, on = c(MarkerName = "rsid"), nomatch = 0
][
    ,.(trait = "ggt", rsid = MarkerName, chr, hg19_pos,
       ea = toupper(ifelse(Effect > 0, Effect_allele, Other_allele)),
       pval = P_value)
]

alt <- fread("../files/alt.txt.gz")

altnoukb <- alt[
    P_value < 5e-8
][
    ref, on = c(MarkerName = "rsid"), nomatch = 0
][
    ,.(trait = "alt", rsid = MarkerName, chr, hg19_pos,
       ea = toupper(ifelse(Effect > 0, Effect_allele, Other_allele)),
       pval = P_value)
]

## AST from Prins et al
## Readme file states column names
ast <- fread("../files/ast.txt.gz",
             col.names = c("rsid", "chr", "hg19_pos",
                           "nea", "ea", "pval", "beta", "se"))

astnoukb <- ast[pval < 5e-8][
   ,.(trait = "ast", rsid, chr, hg19_pos,
      ea = ifelse(beta > 0, ea, nea), pval)
]

## Joining all scores
if(!dir.exists("../files/multiprs_ss/"))
    dir.create("../files/multiprs_ss/")

ls(pattern = "noukb") %>%
    lapply(get) %>%
    lapply(function(x)
        ieugwasr::ld_clump(x, clump_kb = 500, clump_r2 = 0.01,
                           plink_bin = genetics.binaRies::get_plink_binary(),
                           bfile = "../files/1kg_ref/EUR")) %>%
    lapply(function(x) rio::export(x, paste0("../files/multiprs_ss/",
                                             unique(x$trait), ".tsv")))

