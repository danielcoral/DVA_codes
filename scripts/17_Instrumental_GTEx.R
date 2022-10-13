#### GTEx and T2D cross-reference

## Libraries
library(tidyverse)
library(data.table)

## Cross referencing
gtex_bmi_t2d <- fread(paste(
    ## Only SNPs with significant QTL effects in GTEx
    "awk \'FNR==NR{a[$1];next}", ## matrix 'a' with first column (variant_id) as key
    "FNR==1{print $7 FS $3 FS $6;next}", ## preserving column names
    "$1 in a{print $7 FS $3 FS $6}\'", ## values of first column found in matrix a
    "../data/gtex_snps.txt ../data/lookup.txt", "|",
    ## Excluding multi-allelic SNPs
    "awk \'NR==1 || $3==1\'", "|", ## Fourth column is number of alleles per site
    ## Joining with BMI-T2D
    "sed -e '1s/rs_id_dbSNP151_GRCh38p7/rsid/'", "|", ## Changing column name
    "awk 'FNR==NR{a[$1] = $2;next}", ## again, matrix 'a' with first column as key
    "$3 in a{print $0, a[$3]}' OFS=\'\t\' - ../data/bmi_t2d.txt", "|",
    ## Finding SNPs in 1KG reference
    "awk \'NR==FNR{a[$2];next} FNR==1 || ($3 in a)\' ../data/1kg_ref/EUR.bim - "
))[ ## Alignment to BMI increasing allele
   ,`:=`(ea = ifelse(beta.bmi > 0, ea.bmi, nea.bmi),
         nea = ifelse(beta.bmi > 0, nea.bmi, ea.bmi),
         eaf.bmi = ifelse(beta.bmi > 0, eaf.bmi, 1 - eaf.bmi))
][
   ,`:=`(beta.bmi = abs(beta.bmi),
         harmon.t2d = case_when(ea.t2d == ea & nea.t2d == nea ~ 1,
                                ea.t2d == nea & nea.t2d == ea ~ -1, T ~ 0))
][
   ,`:=`(beta.t2d = beta.t2d * harmon.t2d,
         eaf.t2d = ifelse(harmon.t2d == 1, eaf.t2d, 1 - eaf.t2d))
][
    !(harmon.t2d == 0 & abs(eaf.bmi - eaf.t2d) > 0.2)
][
   ,`:=`(eaf = rowMeans(cbind(eaf.bmi, eaf.t2d)),
         pal = "|"(ea %in% c("A", "T") & nea %in% c("A", "T"),
                   ea %in% c("C", "G") & nea %in% c("C", "G")))
][
   ,maf := pmin(eaf, 1 - eaf)
][
    maf > 0.01 & !(pal & maf > .3)
][
   ,.(chr, hg19_pos, hg38_pos = variant_pos, rsid, ea, nea,
      eaf, maf, beta.bmi, se.bmi, p.bmi, n.bmi, beta.t2d, se.t2d, p.t2d, n.t2d)
]

fwrite(gtex_bmi_t2d, "../data/gtex_bmi_t2d.tsv", sep = "\t")

## Update reference panel
dir.create("../data/gtex_bmi_t2d_ref/")

system(
    paste0(
        "cut -d , -f4 ../data/gtex_bmi_t2d.tsv > ../data/gtex_ref_snps.tsv;",
        genetics.binaRies::get_plink_binary(),
        " --bfile ../data/1kg_ref/EUR --extract ../data/gtex_ref_snps.tsv ",
        "--make-bed --out ../data/gtex_bmi_t2d_ref/ref"
    )
)
