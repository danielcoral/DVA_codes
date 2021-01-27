## BMI - T2D instruments

library(tidyverse)

## Intersection between BMI and DIAGRAM (not UKB)
if(!file.exists("../files/bmi_t2diagram.tsv")){
    bmi_t2diagram <- vroom::vroom(
                          ## Using hg19 position
                          pipe(
                              paste0(
                                  "awk \'NR==FNR{a[$1\":\"$2]=$0;next} ",
                                  "$1 in a {print a[$1], \"\t\", $0}\' ",
                                  "../files/bmi.txt ../files/t2d_diagram.txt"
                              )
                          ),
                          col_select = -11,
                          col_names = c("chr", "hg19_pos", "rsid",
                                        "ea.bmi", "nea.bmi", "eaf.bmi",
                                        "beta.bmi", "se.bmi", "p.bmi", "n.bmi",
                                        "x11", # chr:pos (dropping)
                                        "ea.t2d", "nea.t2d",
                                        "beta.t2d", "se.t2d", "p.t2d", "n.t2d")
                      )
    vroom::vroom_write(bmi_t2diagram, "../files/bmi_t2diagram.tsv")
} else {
    bmi_t2diagram <- vroom::vroom("../files/bmi_t2diagram.tsv")
}


## Align to BMI and restrict to variants with significance for BMI
sigbmi_t2d <- bmi_t2diagram %>%
    mutate(
        ## Ensuring allele consistency
        harmon = case_when(
               ea.bmi == ea.t2d & nea.bmi == nea.t2d ~ 1,
               ea.bmi == nea.t2d & nea.bmi == ea.t2d ~ -1,
               T ~ 0
        ),
        ## Aligning to the allele increasing BMI
        ea = ifelse(sign(beta.bmi) == 1, ea.bmi, nea.bmi),
        nea = ifelse(sign(beta.bmi) == 1, nea.bmi, ea.bmi),
        eaf.bmi = ifelse(sign(beta.bmi) == 1, eaf.bmi, 1 - eaf.bmi),
        beta.bmi = abs(beta.bmi),
        beta.t2d = ifelse(ea == ea.t2d, beta.t2d, -beta.t2d),
        ## Removing palindromic SNPs with MAF > 30%
        maf.bmi = ifelse(eaf.bmi <= 0.5, eaf.bmi, 1 - eaf.bmi),
        pal = (ea %in% c("A","T") & nea %in% c("A","T")) |
            (ea %in% c("C", "G") & nea %in% c("C", "G")),
        maf_30 = maf.bmi > 0.3,
        chr_pos = paste(chr, hg19_pos, sep = ":")
    ) %>%
    filter(p.bmi < 5e-8 & !(harmon == 0 | (pal & maf_30) | maf.bmi < 0.01)) %>%
    transmute(chr, pos = hg19_pos, chr_pos, rsid, ea, nea,
              eaf = eaf.bmi,
              beta.bmi, se.bmi, p.bmi, n.bmi,
              beta.t2d, se.t2d, p.t2d, n.t2d,
              disc = ifelse(sign(beta.t2d) == 1, 0, 1))

## Files to write PLINK output
out <- "../files/ld_prox"
targetsname <- paste0(out, ".targets")
outname <- paste0(out, ".targets.ld.gz")
## SNPs to query
utils::write.table(sigbmi_t2d$rsid, file=targetsname,
                   row=FALSE, col=FALSE, qu=FALSE)
## PLINK command
cmd <- paste0(
    genetics.binaRies::get_plink_binary(),
    " --bfile ../files/1kg_ref/EUR", ## Reference
    " --r2 in-phase with-freqs gz",  ## R2 with alleles(PHASE) and MAF
    " --ld-snp-list ", targetsname,  ## SNPs to query
    " --ld-window-kb 1000",          ## 1MB window
    " --ld-window 999999",           ## N of SNPs between each pair
    " --out ", targetsname
)
system(cmd)
## Importing results
ld <- data.table::fread(outname, header = T) %>%
    select(-c(1,2,4)) %>%
    `names<-`(c("orig_rsid", "chr", "pos", "rsid", "phase", "maf", "r2")) %>%
    filter(orig_rsid != rsid) %>%            ## Only proxies
    filter(r2 > 0.5) %>%                     ## r2 threshold
    filter(maf > 0.01) %>%                   ## MAF threshold
    mutate(phase = gsub("/", "", phase)) %>% ## Remove INDELs
    filter(nchar(phase) == 4) %>%
    ## Adding alleles
    {
        df <- .
        alleles_list <- str_extract_all(df$phase, "[A-Z]")
        alleles <- do.call(rbind, alleles_list) %>%
            data.frame() %>%
            `names<-`(c("a1", "b1", "a2", "b2"))
        bind_cols(df, alleles)
    } %>%
    ## Keeping only those not already in the main dataset
    filter(!rsid %in% sigbmi_t2d$rsid)

## Adding association information to proxies
prox_to_use <- left_join(
    mutate(ld, chr_pos = paste(chr, pos, sep = ":")),
    select(sigbmi_t2d, -c(chr, pos, chr_pos)), by = c("orig_rsid" = "rsid")
) %>%
    ## Excluding ambiguous proxies (proxies of both concordant and discordant SNPs)
    {
        df <- .
        prox_to_exclude <- df %>%
            group_by(rsid) %>%
            summarise(discordance = mean(disc)) %>%
            filter(!discordance %in% c(0,1)) %>%
            pull(rsid)
        filter(df, !rsid %in% prox_to_exclude)
    } %>%
    ## Assigning proxies to the highest r2
    group_by(rsid) %>%
    arrange(desc(r2)) %>%
    slice(1) %>%
    ## Aligning to BMI increasing allele
    mutate(ea = ifelse(ea == a1, b1, b2),
           nea = ifelse(nea == a1, b1, b2)) %>%
    select(-c(phase, maf, a1:b2))
        
## Adding proxies to main dataset
sigbmi_t2d <- bind_rows(
    sigbmi_t2d %>%
    mutate(orig_rsid = rsid, r2 = 1) %>%
    select(all_of(names(prox_to_use))),
    prox_to_use
) %>%
    mutate(is_proxy = orig_rsid != rsid)

save(sigbmi_t2d, file = "../files/sigbmi_t2d.RData")
