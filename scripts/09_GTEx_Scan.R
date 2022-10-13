## Processing GTEx data

library(tidyverse)
library(data.table)

qtls <- c("eQTL", "sQTL") ## QTL types available

dir.create("../data/processed_gtex")

gtex_snps <- c() ## Storing SNP ids

## Loop on each dataset
for(q in qtls){
    files <- list.files(paste0("~/projects/DVA/Data/GTEx/GTEx_Analysis_v8_", q),
                        pattern = "signif", full.names = T)
    for(f in files){
        ## Tissue
        tis <- gsub("^.*QTL/|\\.v8.*gz", "", f)
        ## Import each file, MAF > %1, QTL type and tissue indentifiers
        t <- fread(f)[maf > 0.01][,`:=`(qtl_type = q, tissue = tis)][
            ## Extracting chromosome, position, effect and reference alleles
          , c("chr","pos","nea.gtex", "ea.gtex","build") :=
                tstrsplit(variant_id,"_", fixed = T) 
        ][ ## Excluding INDELS
            nchar(ea.gtex) == 1 & nchar(nea.gtex) == 1
        ][ ## Excluding ambiguous palindromic SNPs
            !(((ea.gtex %in% c("A","T") & nea.gtex %in% c("A","T")) & maf > 0.3) |
              ((ea.gtex %in% c("C","G") & nea.gtex %in% c("C","G")) & maf > 0.3))
        ][!chr %in% c("chrX", "chrY")] 
        ## Parsing gene id in sQTL data
        if(any(colnames(t) == "phenotype_id")){
            t[,gene_id := gsub(".*clu.*:", "", phenotype_id)]
            t$phenotype_id <- NULL
        }
        s <- unique(t[,variant_id])
        gtex_snps <- unique(append(gtex_snps, s))
        t %>%
            transmute(tissue, qtl_type, gene_id, chr = as.numeric(gsub("chr", "", chr)),
                      hg38_pos = pos, ea.gtex, nea.gtex, tss_distance, maf, ma_samples, ma_count,
                      slope, slope_se, pval_nominal) %>%
            fwrite(paste0("../data/processed_gtex/", tis, "_", q, ".txt"))
    }
}

### SNP ids of hits
rio::export(data.frame(gtex_snps), "../data/gtex_snps.txt")

## Adding identifiers
system(paste0(
    ## Only SNPs with significant QTL effects in GTEx
    "awk \'NR == FNR{a[$1]; next} FNR==1 || ($1 in a)\' ",
    "../data/gtex_snps.txt ../data/lookup.txt > ../data/gtex_hits.tsv"
))

