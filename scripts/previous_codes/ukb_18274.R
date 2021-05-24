## Additional fields from Naeimeh's application

## UKB application: 18274
origdir <- "/ludc/Active_Projects/UKBB_18274/Private"

newdir <- paste0("/ludc/Home/daniel_c/dva/files/ukb_18274")

if(!dir.exists(newdir))
    dir.create(newdir)

## Additional fields needed
addfields <- list(
    ## For diabetes diagnosis and more fat measurements with DEXA and MRI
    list(basket = 26608,
         fields = c(20008, 20009, 2986, 74, 23245, 23262, 22407, 22408)),
    ## Liver biochemistry
    list(basket = 28686,
         fields = c(30620, 30600, 30650, 30660, 30730, 30840)),
    ## Liver fat
    list(basket3 = 26876,
         fields = 22402)
)

lapply(addfields,
       function(x){
           colsfile <- paste0(newdir, "/cols", x$basket, ".txt")
           enc_file <- paste0(origdir, "/ukb", x$basket, ".enc_ukb")
           write(sort(x$fields), colsfile, ncolumns = 1)
           ## Identifying columns - HTML file
           system(paste("/ludc/Tools/Software/ukbb/ukbconv",
                        enc_file,
                        "docs",
                        "-e/ludc/Tools/Software/ukbb/encoding.ukb",
                        paste0("-i", colsfile),
                        paste0("-o", newdir, "/ukb", x$basket)))
           ## Conversion into R files
           system(paste("/ludc/Tools/Software/ukbb/ukbconv",
                        enc_file,
                        "r",
                        "-e/ludc/Tools/Software/ukbb/encoding.ukb",
                        paste0("-i", colsfile),
                        paste0("-o", newdir, "/ukb", x$basket)))
       })

## Sample file
system(paste("/ludc/Tools/Software/ukbb/ukbgene",
             "imp", "-c1", "-m",
             paste0("-a", paste0(origdir, "/k18274.key"))))

system(paste("mv", list.files(pattern = "\\.sample$"), newdir))

## FAM file
system(paste("/ludc/Tools/Software/ukbb/ukbgene",
             "cal", "-c1", "-m",
             paste0("-a", paste0(origdir, "/k18274.key"))))

system(paste("mv", list.files(pattern = "\\.fam$"), newdir))

system("rm fields.ukb")
