## Insulin DI and sensitivity index (Prokopenko et al)
ins_files <- paste0("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/",
                    c(paste0("ProkopenkoI_24699409_GCST005222/harmonised/",
                             "24699409-GCST005222-EFO_0006832.h.tsv.gz"),
                      paste0("ProkopenkoI_24699409_GCST005221/harmonised/",
                             "24699409-GCST005221-EFO_0004471.h.tsv.gz")))

## Adiponectin (Dastani et al)
adip_url <- paste0("https://www.mcgill.ca/genepi/files/genepi/",
                   "adipogen.discovery.eur_.meta_.public.release.part")
adip_files <- c(paste0(adip_url, 1:5, "_.txt"),
                paste0(adip_url, 6, "__0.txt"))

## Download
destfolder <- "../files/mr_ss/"

download.file(ins_files, method = "libcurl",
              destfile = paste0(destfolder, c("di", "isi"), ".tsv.gz"))

if(!dir.exists(paste0(destfolder, "adip")))
    dir.create(paste0(destfolder, "adip"))

download.file(adip_files, method = "libcurl",
              destfile = paste0(destfolder, "adip/", "adip", 1:6, ".txt"))
