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

## Concordant and discordant SNPs and proxies
system(paste("awk \'NR == FNR{if(FNR == 1){header = $0; next} a[$7] = $0; next}",
             "NR > FNR",
             "{if($2 in a || FNR == 1){printf \"%s\\t%s\\n\", $0, (FNR == 1 ? header : a[$2])}}\'",
             "OFS=\'\t\'",
             "../data/gtex_hits.tsv ../data/proxies.tsv",
             "> ../data/proxies_gtex.tsv"))

proxies_gtex <- fread("../data/proxies_gtex.tsv",
                      drop = c(24,25,26,28,29,30,31,32)) %>%
    rename(hg38_pos = variant_pos)

## Performing the scan
trans_res <- lapply(
    list.files("../data/processed_gtex", full.names = T),
    function(f)
        fread(f)[
          proxies_gtex, on = c("chr", "hg38_pos"), nomatch = 0
        ] %>%
        ## Only 1 association between a SNP and a gene in a tissue
        arrange(ref_rsid, qtl_type, tissue, gene_id, proxy, desc(r2), pval_nominal) %>%
        distinct(ref_rsid, qtl_type, tissue, gene_id, .keep_all = T) %>%
        ## Harmonizing to the BMI increasing allele
        mutate(harmon = case_when(ea.gtex == ea.bmi & nea.gtex == nea.bmi ~ 1,
                                  ea.gtex == nea.bmi & nea.gtex == ea.bmi ~ -1,
                                  TRUE ~ 0),
               slope = slope * harmon)
) %>%
    bind_rows

## Adding gene names
sig_genes <- unique(gsub("\\.[0-9]+$", "", trans_res$gene_id))

ensembl_genes <- biomaRt::useEnsembl(biomart = "ensembl",
                                     dataset = "hsapiens_gene_ensembl")

gene_names <- biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                             filters = "ensembl_gene_id",
                             values = sig_genes,
                             mart = ensembl_genes) %>%
    `names<-`(c("gene_id2", "gene_name"))

## Adding names
trans_res <- trans_res %>%
    mutate(gene_id2 = gsub("\\.[0-9]+", "", gene_id)) %>%
    left_join(gene_names)

fwrite(trans_res, "../data/trans_res.tsv", sep = "\t")

trans_res <- rio::import("../data/trans_res.tsv")
