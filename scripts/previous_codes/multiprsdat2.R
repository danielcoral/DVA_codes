## Preparing PRS

library(googledrive)

if(!dir.exists("../files/multiprs_raw"))
    dir.create("../files/multiprs_raw")

## Data from cohorts other than UK Biobank
## SBP and DBP from the ICBP Consortium
bp.urls <- paste0("http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/",
                  "GCST006001-GCST007000/",
                  c("GCST006259/harmonised/27618452-GCST006259-EFO_0006335.h.tsv.gz",
                    "GCST006258/harmonised/27618452-GCST006258-EFO_0006336.h.tsv.gz"))
                  

download.file(bp.urls[1], "../files/multiprs_raw/sbp_icbp.tsv.gz")

download.file(bp.urls[2], "../files/multiprs_raw/dbp_icbp.tsv.gz")

## HDL, LDL and TG from the GLGC
chol.urls <- paste0("http://csg.sph.umich.edu/abecasis/public/lipids2013/",
                    c("jointGwasMc_LDL.txt.gz", "jointGwasMc_HDL.txt.gz",
                      "jointGwasMc_TG.txt.gz"))

download.file(chol.urls[1], "../files/multiprs_raw/ldl_glgc.txt.gz")
download.file(chol.urls[2], "../files/multiprs_raw/hdl_glgc.txt.gz")
download.file(chol.urls[3], "../files/multiprs_raw/tg_glgc.txt.gz")

## WHR from GIANT
download.file(paste0("https://portals.broadinstitute.org/collaboration/giant/images/5/54/",
                     "GIANT_2015_WHR_COMBINED_EUR.txt.gz"),
              "../files/multiprs_raw/whr_giant.txt.gz")
untar("../files/multiprs_raw/whr_giant.txt.gz", exdir = "../files/multiprs_raw/")
file.rename("../files/multiprs_raw/GIANT_2015_WHR_COMBINED_EUR.txt",
            "../files/multiprs_raw/whr_giant.txt")
unlink("../files/multiprs_raw/whr_giant.txt.gz")

## ALT from Chambers et al.
alt_id <- "1rspqMdM2l143op_OR-I_VUJeFGHCMVTz"
ggt_id <- "1l1ubl7Kj4ryVwG_XJQL5lJGmLvn51YG3"

drive_deauth()
drive_user()
altfile <- drive_get(as_id(alt_id))
ggtfile <- drive_get(as_id(ggt_url))
drive_download(altfile, path = "../files/multiprs_raw/alt.txt.gz")
drive_download(ggtfile, path = "../files/multiprs_raw/ggt.txt.gz")

## AST from Prins et al
options(timeout = 300)
download.file(paste0("http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/",
                     "GCST005001-GCST006000/GCST005064/Prins_28887542_ast.gz"),
              "../files/multiprs_raw/ast.txt.gz")
