## UKB data from project 57232

## Project ID
id <- 57232

## Directory
p_dir <- paste("../files/ukb", id, sep = "_")

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
    ## For QC of samples - genetic sex, PCs, ethnic grouping, kinship...
    22001, 22006, 22021, 22020, 22027, 22019, 20115, 22009, 21000
    ## Age at recruitment, date of assessment, sex, BMI, blood pressure 
    21022, 53, 31, 21001, 4080, 93, 4079, 94,
    ## Vascular or heart problems, non-cancer illnesses, ICD-codes, medications
    6150, 20002, 6154, 6177, 6153, 20003, 41270, 41280, 
    ## Cholesterol, waist and hip circumferences
    30760, 30870, 48, 49, 30690, 30780,
    ## Fat in the extremities measured by impedance
    23120, 23124, 23112, 23116,
    ## Diabetes diagnosis
    2443, 30740, 30750, 2976, 4041
)

write(sort(column_ids), "../files/ukb_57232/column_ids.txt", ncolumns = 1)


## Identifying columns - HTML file
system(paste("/ludc/Tools/Software/ukbb/ukbconv",
             paste(p_dir, paste(ukb_enc, "ukb", sep = "_"), sep = "/"),
             "docs",
             "-e/ludc/Tools/Software/ukbb/encoding.ukb",
             "-i../files/ukb_57232/column_ids.txt"))

## Conversion into R files
system(paste("/ludc/Tools/Software/ukbb/ukbconv",
             paste(p_dir, paste(ukb_enc, "ukb", sep = "_"), sep = "/"),
             "r",
             "-e/ludc/Tools/Software/ukbb/encoding.ukb",
             "-i../files/ukb_57232/column_ids.txt"))

## Sample file
system(paste("/ludc/Tools/Software/ukbb/ukbgene",
             "imp", "-c1", "-m",
             paste0("-a", paste(p_dir, ukb_key, sep = "/"))))

system(paste("mv", list.files(pattern = "\\.sample$"), p_dir))

## FAM file
system(paste("/ludc/Tools/Software/ukbb/ukbgene",
             "cal", "-c5", "-m",
             paste0("-a", paste(p_dir, ukb_key, sep = "/"))))

system(paste("mv", list.files(pattern = "\\.fam$"), p_dir))

## Moving fields file to directory
system(paste("mv fields.ukb", p_dir))
