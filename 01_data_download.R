#### Data download

### Creating directories needed
if(!dir.exists("../files"))
    dir.create("../files")

if(!dir.exists("../plots"))
    dir.create("../plots")


### Downloading data needed

## BMI from the GIANT Consortium + UK Biobank
if(!file.exists("../files/bmi.txt")){
    download.file(
        paste0(
            "https://portals.broadinstitute.org/collaboration/giant/images/c/c8/",
            "Meta-analysis_Locke_et_al%2BUKBiobank_2018_UPDATED.txt.gz"
        ),
        "../files/bmi.txt.gz"
    )
    R.utils::gunzip("../files/bmi.txt.gz")
}

## T2D from DIAGRAM Consortium + UK Biobank
if(!file.exists("../files/t2d.txt")){
    system(
        paste0(
            'curl --data "confirm=1" --output ',
            '../files/t2d.zip',
            ' https://www.diagram-consortium.org/check_check.Mahajan2018b.T2D.3.php'
        )
    )
    unzip("../files/t2d.zip", exdir = "../files/")
    unlink("../files/t2d.zip")
    file.rename(
        "../files/Mahajan.NatGenet2018b.T2D.European.txt",
        "../files/t2d.txt"
    )
}

## T2D from DIAGRAM only - for MR
if(!file.exists("../files/t2d_diagram.txt")){
    system(
        paste0(
            'curl --data "confirm=1" --output ',
            '../files/t2d_diagram.zip',
            ' https://www.diagram-consortium.org/check_check.scott_1.php'
        )
    )
    unzip("../files/t2d_diagram.zip", exdir = "../files/")
    unlink("../files/t2d_diagram.zip")
    file.rename(
        "../files/METAANALYSIS_DIAGRAM_SE1.txt",
        "../files/t2d_diagram.txt"
    )
}    

## Reference panel for clumpling: 1KG without palindromic or INDELs
if(!dir.exists("../files/1kg_ref"))
    dir.create("../files/1kg_ref")

if(!file.exists("../files/1kg.tgz"))
    download.file("http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz",
                  "../files/1kg.tgz")

if(!identical(list.files("../files/1kg_ref"),
              c("EUR.bed", "EUR.bim", "EUR.fam")))
    untar("../files/1kg.tgz",
          files = c("EUR.bed", "EUR.bim", "EUR.fam"),
          exdir = "../files/1kg_ref/")

