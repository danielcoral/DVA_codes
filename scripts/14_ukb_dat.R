## UKB data from project 57232

## Project ID
id <- 57232

## Directory
p_dir <- paste("/ludc/Home/daniel_c/dva/files/ukb", id, sep = "_")

if(!dir.exists(p_dir))
    dir.create(p_dir)

## Files needed in this directory:
## 1. The key (sent through e-mail)
## 2. The .enc file (manual download)

## This project only has one basket (41806)

basket <- 41806

ukb_key <- paste0("k", id, "r", basket, ".key")

ukb_enc <- paste0("ukb", basket, ".enc")

## Making sure the dataset is correctly downloaded
system(paste("/ludc/Tools/Software/ukbb/ukbmd5",
             paste(p_dir, ukb_enc, sep = "/")))

## Decompressing files
## This needs the helper programs already in the server

## We need the keyvalue inside the key file - second line
keyvalue <- readLines(paste(p_dir, ukb_key, sep = "/"))[2]

system(paste("/ludc/Tools/Software/ukbb/ukbunpack",
             paste(p_dir, ukb_enc, sep = "/"), keyvalue))

## Columns needed, based on clusters found
## Selected through the Data Showcase
column_ids <- c(
    ## For QC of samples:
    ## Genetic sex, ethnic grouping, used in genetic PC (unrelated samples),
    ## Heterozygosity outliers, sex chromosome aneuploidy, genetic PCs
    22001, 22006, 22020, 22027, 22019, 22009,
    ## Age at recruitment, date of assessment, sex
    21022, 53, 31
)

write(sort(column_ids), "~/dva/files/ukb_57232/column_ids.txt", ncolumns = 1)

## Identifying columns - HTML file
system(paste("/ludc/Tools/Software/ukbb/ukbconv",
             paste(p_dir, paste(ukb_enc, "ukb", sep = "_"), sep = "/"),
             "docs",
             "-e/ludc/Tools/Software/ukbb/encoding.ukb",
             "-i/ludc/Home/daniel_c/dva/files/ukb_57232/column_ids.txt"))

## Conversion into R files
system(paste("/ludc/Tools/Software/ukbb/ukbconv",
             paste(p_dir, paste(ukb_enc, "ukb", sep = "_"), sep = "/"),
             "r",
             "-e/ludc/Tools/Software/ukbb/encoding.ukb",
             "-i/ludc/Home/daniel_c/dva/files/ukb_57232/column_ids.txt"))

## Sample file
system(paste("/ludc/Tools/Software/ukbb/ukbgene",
             "imp", "-c1", "-m",
             paste0("-a", paste(p_dir, ukb_key, sep = "/"))))

system(paste("mv", list.files(pattern = "\\.sample$"), p_dir))

## FAM file
system(paste("/ludc/Tools/Software/ukbb/ukbgene",
             "cal", "-c1", "-m",
             paste0("-a", paste(p_dir, ukb_key, sep = "/"))))

system(paste("mv", list.files(pattern = "\\.fam$"), p_dir))

## Moving fields file to directory
system(paste("mv fields.ukb", p_dir))
