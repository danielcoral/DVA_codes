## Data from MRC IEU database
ids_ieu <- c(
    ## The betas of all these variables are in SD units
    armfatmass_l = "ukb-b-8338", armfatmass_r = "ukb-b-6704",
    legfatmass_l = "ukb-b-7212", legfatmass_r = "ukb-b-18096",
    hdl = "ukb-d-30760_irnt", ldl = "ukb-d-30780_irnt",
    tg = "ukb-d-30870_irnt"
)

options(timeout = 300)
lapply(ids_ieu,
       function(x){
           download.file(
               paste0("https://gwas.mrcieu.ac.uk/files/",
                      paste0(x, "/", x, ".vcf.gz")),
               destfile = paste0("../files/mr_ss/", x, ".vcf.gz")
           )
           download.file(
               paste0("https://gwas.mrcieu.ac.uk/files/",
                      paste0(x, "/", x, ".vcf.gz.tbi")),
               destfile = paste0("../files/mr_ss/", x, ".vcf.gz.tbi")
           )
       })

## Data for WHR from GIANT + UKB
whr_url <- paste0("https://zenodo.org/record/1251813/files/",
                  "whr.giant-ukbb.meta-analysis.combined.23May2018.txt.gz?download=1")
download.file(whr_url, destfile = "../files/mr_ss/whr.txt.gz")

## From GWAS Catalog
gwascat_ftp <- "ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"

gwascat_urls <- c(
    ## Blood pressure (in natural units)
    sbp = "EvangelouE_30224653_GCST006624/Evangelou_30224653_SBP.txt.gz",
    dbp = "EvangelouE_30224653_GCST006630/Evangelou_30224653_DBP.txt.gz"
)

lapply(names(gwascat_urls),
       function(x) download.file(paste0(gwascat_ftp, gwascat_urls[[x]]),
                                 destfile = paste0("../files/mr_ss/", x, ".txt.gz")))
