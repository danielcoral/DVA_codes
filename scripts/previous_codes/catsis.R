## Idea de catsis

library(tidyverse)

bmi <- data.table::fread("../files/bmi.txt")

bmisig <- bmi %>%
    rename(eaf = Freq_Tested_Allele_in_HRS,
           ea = Tested_Allele, nea = Other_Allele,
           rsid = SNP, pval = P) %>%
    filter(pval < 5e-8, eaf > 0.01, eaf < 0.99) %>%
    mutate(pal1 = ea %in% c("A", "T") & nea %in% c("A", "T"),
           pal2 = ea %in% c("C", "G") & nea %in% c("C", "G"),
           maf = ifelse(eaf > 0.5, 1 - eaf, eaf),
           ambig = (pal1 | pal2) & maf > 0.3) %>%
    ieugwasr::ld_clump(clump_kb = 500, clump_r2 = 0.01,
                       bfile = "../files/1kg_ref/EUR",
                       plink_bin = genetics.binaRies::get_plink_binary()) %>%
    mutate(ea2 = ifelse(BETA > 0, ea, nea),
           nea2 = ifelse(BETA > 0, nea, ea),
           ea = ea2, nea = nea2,
           beta = abs(BETA),
           hg19_pos = POS) %>%
    select(-c(POS, ea2, nea2, BETA, pal1:id)) %>%
    rename_with(tolower)
    

bmisig %>%
    head
   
