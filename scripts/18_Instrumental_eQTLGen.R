## Blood eQTL MR

library(tidyverse)
library(ieugwasr)
library(ggh4x)
## SMR & HEIDI functions
source(paste0("./", list.files(pattern = "MR_FX.R")))

## Blood eQTL hits
bloodeqtl_scan <- rio::import("../data/bloodeqtl_scan.tsv")

## Genes
sig_genes <- unique(bloodeqtl_scan$id)

## Instruments
bloodeqtl_ins <- tophits(sig_genes, clump = 0)

rio::export(bloodeqtl_ins, "../data/bloodeqtl_ins.tsv")

bloodeqtl_ins <- rio::import("../data/bloodeqtl_ins.tsv")

## BMI and T2D data
bmi_t2d <- rio::import("../data/bmi_t2d.txt") %>%
    mutate(ea = ifelse(beta.bmi > 0, ea.bmi, nea.bmi),
           nea = ifelse(beta.bmi > 0, nea.bmi, ea.bmi),
           eaf.bmi = ifelse(beta.bmi > 0, eaf.bmi, 1 - eaf.bmi),
           beta.bmi = abs(beta.bmi),
           harmon = case_when(ea.t2d == ea & nea.t2d == nea ~ 1,
                              nea.t2d == ea & ea.t2d == nea ~ -1),
           beta.t2d = beta.t2d * harmon,
           eaf.t2d = ifelse(harmon == 1, eaf.t2d, 1 - eaf.t2d)) %>%
    filter(!is.na(harmon), abs(eaf.bmi - eaf.t2d) < 0.2) %>%
    mutate(eaf = rowMeans(cbind(eaf.bmi, eaf.t2d)),
           maf = ifelse(eaf <= .5, eaf, 1 - eaf),
           ambig = "&"((ea %in% c("A", "T") & nea %in% c("A", "T") |
                        ea %in% c("C", "G") & nea %in% c("C", "G")),
                       maf > .3)) %>%
    filter(!ambig) %>%
    select(-c(harmon, eaf.bmi, eaf.t2d, ambig, ea.bmi, nea.bmi, ea.t2d, nea.t2d))

## Cis regions
gene_ids <- unique(bloodeqtl_ins$trait)

## To look for hg19 genomic locations
ensembl_genes_grch37 <- biomaRt::useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL",
                                            dataset = "hsapiens_gene_ensembl",
                                            host = "grch37.ensembl.org")

## SMR
bloodeqtl_smr <- biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name",
                                        "chromosome_name",
                                        "start_position", "end_position",
                                        "transcript_start"),
                         filters = "ensembl_gene_id",
                         values = gene_ids, mart = ensembl_genes_grch37) %>%
    group_by(ensembl_gene_id) %>%
    summarise(gene_name = unique(external_gene_name),
              chr = as.numeric(unique(chromosome_name)),
              min_pos = min(transcript_start - 500000, start_position, end_position),
              max_pos = max(transcript_start + 500000, start_position, end_position)) %>%
    filter(!is.na(chr)) %>%
    rename(trait = ensembl_gene_id) %>%
    inner_join(mutate(bloodeqtl_ins, chr = as.numeric(chr))) %>%
    filter(position > min_pos, position < max_pos) %>%
    rename(hg19_pos = position) %>%
    rename_with(~paste(.x, "exp", sep = "."),
                c(ea, nea, eaf, beta, se, p, n)) %>%
    inner_join(bmi_t2d) %>%
    filter(p.bmi < 5e-8) %>%
    mutate(harmon = case_when(ea == ea.exp & nea == nea.exp ~ 1,
                              ea == nea.exp & nea == ea.exp ~ -1),
           beta.exp = beta.exp * harmon, eaf.exp = ifelse(harmon == 1, eaf.exp, 1 - eaf.exp)) %>%
    filter(!is.na(harmon), abs(eaf - eaf.exp) < .2) %>%
    select(-c(harmon, ea.exp, nea.exp, eaf.exp)) %>%
    group_by(trait) %>%
    arrange(p.exp) %>% slice(1) %>%
    ungroup %>%
    mutate(beta.smr = beta.t2d/beta.exp,
           se.smr = std_err_smr(beta.exp, se.exp, beta.t2d, se.t2d),
           p.smr = pval_smr(beta.exp, se.exp, beta.t2d, se.t2d),
           padj = p.adjust(p.smr, "fdr"),
           disc = sign(beta.t2d/beta.bmi))

rio::export(bloodeqtl_smr, "../data/bloodeqtl_smr.tsv")

bloodeqtl_smr <- rio::import("../data/bloodeqtl_smr.tsv")

bloodeqtl_sigsmr <- filter(bloodeqtl_smr, padj < 0.05)

bloodeqtl_loci <- bloodeqtl_sigsmr %>%
    group_by(id) %>% 
    group_modify(~{
        cisreg <- filter(bmi_t2d, chr == .x$chr,
                         hg19_pos > (.x$hg19_pos - 500000),
                         hg19_pos < (.x$hg19_pos + 500000))
        associations(cisreg$rsid, id = .y$id, proxies = 0) %>%
            select(-c(id, chr, position)) %>%
            rename_with(~paste(.x, "exp", sep = "."),
                        c(beta, se, p, n, ea, nea, eaf)) %>%
            inner_join(cisreg) %>%
            mutate(top = rsid == .x$rsid,
                   harmon = case_when(ea.exp == ea & nea.exp == nea ~ 1,
                                      ea.exp == nea & nea.exp == ea ~ -1),
                   beta.exp = beta.exp * harmon,
                   eaf.exp = ifelse(harmon == 1, eaf.exp, 1 - eaf.exp)) %>%
            filter(!is.na(harmon), abs(eaf - eaf.exp) < .2) %>%
            select(-c(harmon, ea.exp, nea.exp, ea.exp))
    }) %>%
    ungroup

rio::export(bloodeqtl_loci, "../data/bloodeqtl_loci.tsv")

bloodeqtl_loci <- rio::import("../data/bloodeqtl_loci.tsv")

## Making reference panel for clumping
plink <- genetics.binaRies::get_plink_binary()
bfiledir <- "/ludc/Home/daniel_c/dva/files/bloodeqtl_ins_ref"
dir.create(bfiledir)
rio::export(unique(select(bloodeqtl_loci, rsid)), "../data/bloodeqtl_snps.tsv")
system(paste(plink,
             "--bfile ../data/1kg_ref/EUR",
             "--extract ../data/bloodeqtl_snps.tsv",
             "--make-bed --out", paste(bfiledir, "ref", sep = "/")))

## HEIDI
bloodeqtl_heidi <- bloodeqtl_loci %>%
    group_by(id) %>% 
    group_modify(~{
        ## To be returned if HEIDI fails
        nullres <- data.frame(p.heidi.bmi = NaN, p.heidi.t2d = NaN)
        ## Enough SNPs for HEIDI
        if(nrow(.x) >= 3) {
            ## LD matrix
            ldmat <- tryCatch(ld_matrix(.x$rsid, bfile = paste(bfiledir, "ref", sep = "/"),
                                        plink_bin = plink),
                              error = function(e) NULL)
            no_ldmat <- is.null(ldmat) | is.null(dim(ldmat))
            if(!no_ldmat) {
                ## Correcting error with Thymine (read as TRUE)
                if(any(grepl("TRUE", rownames(ldmat))))
                    colnames(ldmat) <- rownames(ldmat) <- gsub("TRUE","T", rownames(ldmat))
                ## Choosing instruments for HEIDI
                ldmat_snps <- tibble(rsid = rownames(ldmat)) %>%
                    separate(rsid, c("rsid", "a2", "a1"), sep = "_") %>%
                    rowid_to_column()
                heidi_ins <- inner_join(.x, ldmat_snps, by = "rsid") %>%
                    mutate(harmon = case_when(ea == a1 & nea == a2 ~ 1,
                                              ea == a2 & nea == a1 ~ -1,
                                              T ~ 0)) %>%
                    filter(harmon != 0)
                ## In case instruments could not be harmonized
                heidi_ins <- tryCatch(
                    mutate(heidi_ins,
                           harmon = ifelse(harmon[top] == 1, harmon, -harmon),
                           r2 = ldmat[rowid[top], rowid]^2) %>%
                    ## Top 20 instruments with r2 < 0.9 and > 0.05
                    filter(top | (r2 > 0.05 & r2 < 0.9)) %>%
                    arrange(desc(top), desc(r2), p.exp) %>% slice(1:20),
                    error = function(e) data.frame(notrun = 1)
                )
                ## Minimum number of instruments needed for HEIDI = 3
                if(nrow(heidi_ins) >= 3){
                    ## Align LD matrix
                    ldmat <- ldmat[heidi_ins$rowid, heidi_ins$rowid]
                    for(i in 1:nrow(ldmat)){
                        for(j in 1:nrow(ldmat)){
                            ldmat[i,j] <- ifelse(identical(heidi_ins$harmon[i],
                                                           heidi_ins$harmon[j]),
                                                 ldmat[i,j], -ldmat[i,j])
                        }
                    }
                    ## Calculating HEIDI p-value
                    data.frame(p.heidi.bmi = tryCatch(heidi_pvalue(heidi_ins$beta.exp,
                                                                   heidi_ins$se.exp,
                                                                   heidi_ins$beta.bmi,
                                                                   heidi_ins$se.bmi,
                                                                   ldmat, 1),
                                                      error = function(e) NaN),
                               p.heidi.t2d = tryCatch(heidi_pvalue(heidi_ins$beta.exp,
                                                                   heidi_ins$se.exp,
                                                                   heidi_ins$beta.t2d,
                                                                   heidi_ins$se.t2d,
                                                                   ldmat, 1),
                                                      error = function(e) NaN))
                } else { message("< 3 SNPs for HEIDI"); return(nullres) }
            } else { message("Unvalid LD matrix for probe"); return(nullres) }
        } else { message("< 3 SNP in found in the cis region"); return(nullres) }
    }) %>%
    ungroup %>%
    filter(!is.na(p.heidi.bmi), !is.na(p.heidi.t2d))

bloodeqtl_smrheidi <- inner_join(bloodeqtl_sigsmr, bloodeqtl_heidi)

rio::export(bloodeqtl_smrheidi, "../data/bloodeqtl_smrheidi.tsv")

