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

## SMR functions
source(paste0("./", list.files(pattern = "MR_FX.R")))

### Transcriptome scan
trans_res <- fread("../data/trans_res.tsv")

## GTEx-BMI-T2D cross-reference
gtex_bmi_t2d <- fread("../data/gtex_bmi_t2d.tsv") %>%
    rename(maf_bmit2d = maf)

bfiledir <- "/ludc/Home/daniel_c/dva/files/gtex_bmi_t2d_ref/ref"

## Implementing SMR
smr_res <- trans_res %>%
    group_by(tissue, qtl_type) %>%
    group_modify(~{
        filenam <- paste0("../data/processed_gtex/",
                          paste(.y$tissue, .y$qtl_type, sep = "_"), ".txt")
        fread(filenam)[
            gene_id %in% unique(.x$gene_id) & pval_nominal < 5e-8
        ][
            gtex_bmi_t2d, on = c("chr", "hg38_pos"), nomatch = 0
        ][
            p.bmi < 5e-8 & abs(maf - maf_bmit2d) < 0.2
        ][
           ,`:=`(harmon = case_when(ea == ea.gtex & nea == nea.gtex ~ 1,
                                    ea == nea.gtex & nea == ea.gtex ~ 1, T ~ 0))
        ][
            harmon != 0
        ][
           ,slope := slope * harmon
        ][
           ,retain := if(.N > 1) {
                          plink <- genetics.binaRies::get_plink_binary()
                          dat <- tibble(rsid, pval = pval_nominal)
                          cdat <- ieugwasr::ld_clump(dat, clump_kb = 500, clump_r2 = 0.01,
                                                     bfile = bfiledir,
                                                     plink_bin = plink)
                          rsid %in% cdat$rsid
                      } else TRUE,
            by = gene_id
        ][
            retain == TRUE
        ][
           ,c("retain", "tissue", "qtl_type", "harmon", "ea.gtex", "nea.gtex") := NULL
        ]
    }) %>%
    ungroup %>%
    unique %>%
    mutate(beta.smr = beta.t2d/slope,
           se.smr = std_err_smr(slope, slope_se, beta.t2d, se.t2d),
           p.smr = pval_smr(slope, slope_se, beta.t2d, se.t2d),
           disc = sign(beta.t2d/beta.bmi),
           padj = p.adjust(p.smr, "fdr"))

rio::export(smr_res, "../data/smr_res.tsv")

## Significant SMR results
smr_sig <- fread("../data/smr_res.tsv") %>%
    filter(padj < 0.05)

## GTEx-BMI-T2D cross-reference
gtex_bmi_t2d <- fread("../data/gtex_bmi_t2d.tsv") %>%
    rename(maf_bmi_t2d = maf)

## Implementation
heidi_res <- smr_sig %>%
    group_by(qtl_type, tissue) %>%
    group_modify(~{
        ## 1. QTL data for each tissue
        filenam <- paste0("../data/processed_gtex/",
                          paste(.y$tissue, .y$qtl_type, sep = "_"), ".txt")
        dat <- fread(filenam)[
            ## 2. Extracting information for the genes of interest
            gene_id %in% unique(.x$gene_id)
        ][
            ## 3. Adding BMI and T2D data
            gtex_bmi_t2d, on = c("chr", "hg38_pos"), nomatch = 0
        ][
            ## 4. Harmonization
            abs(maf - maf_bmi_t2d) < 0.2
        ][
           ,`:=`(harmon = case_when(ea == ea.gtex & nea == nea.gtex ~ 1,
                                    ea == nea.gtex & nea == ea.gtex ~ 1, T ~ 0))
        ][
            harmon != 0
        ][
           ,`:=`(slope = slope * harmon)
        ] %>% unique
        ## 5. Extracting information relevant for each probe
        lapply(1:nrow(.x),
               function(i){
                   probe <- data.frame(.x[i,])
                   ## 6. Instruments in the cis region of the probe (at least 3 needed)
                   cis_reg <- filter(dat, gene_id == probe$gene_id,
                                     hg19_pos > probe$hg19_pos - 500000,
                                     hg19_pos < probe$hg19_pos + 500000)
                   if(nrow(cis_reg) >= 3) {
                       ## 7. LD matrix
                       plink <- genetics.binaRies::get_plink_binary()
                       b.file <- "/ludc/Home/daniel_c/dva/files/gtex_bmi_t2d_ref/ref"
                       ldmat <- tryCatch(ieugwasr::ld_matrix(cis_reg$rsid, bfile = b.file,
                                                             plink_bin = plink),
                                         error = function(e) NULL)
                       no_ldmat <- is.null(ldmat) | is.null(dim(ldmat))
                       if(!no_ldmat) {
                           ## 8. Correcting error with Thymine (read as TRUE)
                           if(any(grepl("TRUE", rownames(ldmat))))
                               colnames(ldmat) <- rownames(ldmat) <- gsub("TRUE","T",
                                                                          rownames(ldmat))
                           ## 9. Choosing instruments for HEIDI
                           ldmat_snps <- tibble(rsid = rownames(ldmat)) %>%
                               separate(rsid, c("rsid", "a2", "a1"), sep = "_") %>%
                               rowid_to_column()
                           heidi_ins <- inner_join(cis_reg, ldmat_snps, by = c("rsid")) %>%
                               mutate(harmon = case_when(ea == a1 & nea == a2 ~ 1,
                                                         ea == a2 & nea == a1 ~ -1,
                                                         T ~ 0)) %>% filter(harmon != 0) %>%
                               mutate(top = rsid == probe$rsid,
                                      harmon = ifelse(harmon[top] == 1, harmon, -harmon),
                                      r2 = ldmat[rowid[top], rowid]^2) %>%
                               ## Top 20 instruments with r2 < 0.9 and > 0.05
                               filter(top | (r2 > 0.05 & r2 < 0.9)) %>%
                               arrange(desc(top), desc(r2), pval_nominal) %>% slice(1:20)
                           ## 9. Minimum number of instruments needed for HEIDI = 3
                           if(nrow(heidi_ins) >= 3){
                               ## 10. Align LD matrix
                               ldmat <- ldmat[heidi_ins$rowid, heidi_ins$rowid]
                               for(i in 1:nrow(ldmat)){
                                   for(j in 1:nrow(ldmat)){
                                       ldmat[i,j] <- ifelse(identical(heidi_ins$harmon[i],
                                                                      heidi_ins$harmon[j]),
                                                            ldmat[i,j], -ldmat[i,j])
                                   }
                               }
                               ## 11. Calculating HEIDI p-value for the probe
                               mutate(probe,
                                      p.heidi.bmi = tryCatch(heidi_pvalue(heidi_ins$slope,
                                                                          heidi_ins$slope_se,
                                                                          heidi_ins$beta.bmi,
                                                                          heidi_ins$se.bmi,
                                                                          ldmat, 1),
                                                             error = function(e) NULL),
                                      p.heidi.t2d = tryCatch(heidi_pvalue(heidi_ins$slope,
                                                                          heidi_ins$slope_se,
                                                                          heidi_ins$beta.t2d,
                                                                          heidi_ins$se.t2d,
                                                                          ldmat, 1),
                                                             error = function(e) NULL))
                           } else { message(paste("< 3 SNPs for HEIDI for probe", probe$rsid,
                                                  "in gene", probe$gene_id, "-",
                                                  .y$tissue, "-", .y$qtl_type)); return(NULL) }
                       } else { message(paste("Unvalid LD matrix for probe", probe$rsid,
                                              "in gene", probe$gene_id, "-",
                                              .y$tissue, "-", .y$qtl_type)); return(NULL) }
                   } else { message(paste("< 3 SNP in found in the cis region for probe",
                                          probe$rsid, "in gene", probe$gene_id, "-",
                                          .y$tissue, "-", .y$qtl_type)); return(NULL) }
               }) %>% bind_rows
    })

## Saving
rio::export(heidi_res, "../data/heidi_res.tsv")

