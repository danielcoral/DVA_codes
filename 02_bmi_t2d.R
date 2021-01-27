#### Intersection and discordance

### The tidyverse
library(tidyverse)
library(data.table)
library(biomaRt)

### BMI - T2D Cross-referencing
if(!file.exists("../files/bmi_t2d.txt")){
    bmi_t2d <- vroom::vroom(
                          ## Using hg19 position
                          pipe(
                              paste0(
                                  "awk \'NR==FNR{a[$1\":\"$2]=$0;next} ",
                                  "$1 in a {print a[$1], \"\t\", $0}\' ",
                                  "../files/bmi.txt ../files/t2d.txt"
                              )
                          ),
                          col_select = -c(11,12,13),
                          col_names = c("chr", "hg19_pos", "rsid",
                                        "ea.bmi", "nea.bmi", "eaf.bmi",
                                        "beta.bmi", "se.bmi", "p.bmi", "n.bmi",
                                        ## Variables to drop
                                        "x1", "x2", "x3",
                                        "ea.t2d", "nea.t2d", "eaf.t2d",
                                        "beta.t2d", "se.t2d", "p.t2d", "n.t2d")
                      )
    vroom::vroom_write(bmi_t2d, "../files/bmi_t2d.txt")
    setDT(bmi_t2d)
} else {
    bmi_t2d <- vroom::vroom("../files/bmi_t2d.txt")
    setDT(bmi_t2d)
}


## Genome-wide significant for both conditions
mix <- bmi_t2d[
    p.bmi <= 5e-8 & p.t2d <= 5e-8
][,
  ## Allele consistency
  harmon := fcase(
      ea.bmi == ea.t2d & nea.bmi == nea.t2d, 1,
      ea.bmi == nea.t2d & nea.bmi == ea.t2d, -1,
      default = 0
  )
  ][harmon != 0
    ][,
      ## Aligning to the BMI increasing allele
      `:=`(
          ea = ifelse(sign(beta.bmi) == 1, ea.bmi, nea.bmi),
          nea = ifelse(sign(beta.bmi) == 1, nea.bmi, ea.bmi),
          eaf.bmi = ifelse(sign(beta.bmi) == 1, eaf.bmi, 1 - eaf.bmi)
          )
      ][,
        `:=`(
            beta.bmi = abs(beta.bmi),
            beta.t2d = ifelse(ea == ea.t2d, beta.t2d, -beta.t2d),
            eaf.t2d = ifelse(ea == ea.t2d, eaf.t2d, 1 - eaf.t2d),
            maf.bmi = ifelse(eaf.bmi <= 0.5, eaf.bmi, 1 - eaf.bmi)
        )
        ][,
          maf.t2d := ifelse(eaf.t2d <= 0.5, eaf.t2d, 1 - eaf.t2d)
          ][,
            ## Removing rare SNPs & unequal MAF & palindromic SNPs with MAF > 0.3
            `:=`(
                rare = maf.bmi < 0.01 | maf.t2d < 0.01,
                uneq_maf = abs(maf.bmi - maf.t2d) > 0.2,
                pal = (ea %in% c("A","T") & nea %in% c("A","T")) |
                    (ea %in% c("C", "G") & nea %in% c("C", "G")),
                maf_30 = maf.bmi > 0.3 | maf.t2d > 0.3
            )
            ][!(pal & maf_30 | rare | uneq_maf)] %>%
    { ## Clumping based on BMI p values
        dat <- transmute(., rsid, pval = p.bmi)
        clump_dat <- ieugwasr::ld_clump(
                                   dat,
                                   clump_kb = 500,
                                   clump_r2 = 0.01,
                                   plink_bin = genetics.binaRies::get_plink_binary(),
                                   bfile = "../files/1kg_ref/EUR"
                                   )
        .[rsid %in% clump_dat$rsid] } %>%
    transmute(chr, hg19_pos, rsid, ea, nea,
              eaf.bmi, beta.bmi, se.bmi, p.bmi, n.bmi,
              eaf.t2d, beta.t2d, se.t2d, p.t2d, n.t2d,
              disc = ifelse(sign(beta.t2d) == 1, 0, 1))

## Annotating SNPs to genes using Ensembl Biomart
ensembl_snps <- useEnsembl(biomart = "snps", dataset = "hsapiens_snp")
ensembl_genes <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

nearest_genes <- lapply(
    mix$rsid,
    function(x){
        ## Gene containing the SNP
        snp_info <- getBM(attributes = c("refsnp_id", "chr_name", "chrom_start", "chrom_end",
                                         "ensembl_gene_stable_id", "ensembl_type"),
                          filters = "snp_filter",
                          values = x,
                          mart = ensembl_snps) %>%
            ## Removing duplicates
            distinct(refsnp_id, .keep_all = T) %>%
            ## Vicinity of the SNP (1MB window) to look for genes
            mutate(chr_reg = paste(chr_name, chrom_start - 1000000, chrom_end + 1000000, sep = ":")) 
        if(is.na(snp_info$ensembl_gene_stable_id)){
            close_genes <- getBM(attributes = c("external_gene_name", "chromosome_name", "strand",
                                                "start_position", "end_position", "gene_biotype"),
                                 filters = "chromosomal_region",
                                 values = snp_info$chr_reg,
                                 mart = ensembl_genes) %>%
                mutate(rsid = x,
                       snp_pos = snp_info$chrom_start,
                       ## Distance from the beginning of the gene
                       tss = ifelse(strand == 1, start_position, end_position),
                       tss_dist = abs(snp_pos - tss),
                       ## Criteria to prioritize genes to assign to the SNP:
                       ## 1. Protein-coding gene
                       protcod = gene_biotype == "protein_coding",
                       ## 2. Antisense of protein-coding
                       as_gene = grepl("-AS", external_gene_name),
                       ## 3. Variant inside a long intergenic non-coding
                       inside_gene = snp_pos > start_position & snp_pos < end_position,
                       inside_linc = inside_gene & grepl("linc", external_gene_name, ignore.case = T),
                       ## These three are treated as equally important
                       priority = protcod | as_gene | inside_linc) %>%
                arrange(desc(priority), desc(inside_gene), tss_dist) %>%
                slice(1)
            res <- close_genes
        } else {
            closest_gene <- getBM(attributes = "external_gene_name",
                                  filters = "ensembl_gene_id",
                                  values = snp_info$ensembl_gene_stable_id,
                                  mart = ensembl_genes) %>%
                mutate(rsid = x)
            res <- closest_gene
        }
        message(paste(x, "done"))
        return(dplyr::select(res, rsid, nearest_gene = external_gene_name))
    }
) %>%
    bind_rows

head(nearest_genes)

mix <- left_join(mix, nearest_genes, by = "rsid")

vroom::vroom_write(mix, "../files/mix.txt")
