#### Data download

### Creating directories needed
if(!dir.exists("~/dva/files"))
    dir.create("~/dva/files")

### Downloading data needed

## BMI from the GIANT Consortium + UK Biobank
if(!file.exists("~/dva/files/bmi.txt")){
    download.file(
        paste0(
            "https://portals.broadinstitute.org/collaboration/giant/images/c/c8/",
            "Meta-analysis_Locke_et_al%2BUKBiobank_2018_UPDATED.txt.gz"
        ),
        "~/dva/files/bmi.txt.gz"
    )
    R.utils::gunzip("~/dva/files/bmi.txt.gz")
}

## T2D from DIAGRAM Consortium + UK Biobank
if(!file.exists("~/dva/files/t2d.txt")){
    system(
        paste0(
            'curl --data "confirm=1" --output ',
            '~/dva/files/t2d.zip',
            ' https://www.diagram-consortium.org/check_check.Mahajan2018b.T2D.3.php'
        )
    )
    unzip("~/dva/files/t2d.zip", exdir = "~/dva/files/")
    unlink("~/dva/files/t2d.zip")
    file.rename(
        "~/dva/files/Mahajan.NatGenet2018b.T2D.European.txt",
        "~/dva/files/t2d.txt"
    )
}

## BMI from GIANT only
if(!file.exists("~/dva/files/bmi_giant.txt")){
    download.file(
        paste0(
            "https://portals.broadinstitute.org/collaboration/giant/images/1/15/",
            "SNP_gwas_mc_merge_nogc.tbl.uniq.gz"
        ),
        "~/dva/files/bmi_giant.txt.gz"
    )
    R.utils::gunzip("~/dva/files/bmi_giant.txt.gz")
}

## T2D from DIAGRAM only
if(!file.exists("~/dva/files/t2d_diagram.txt")){
    system(
        paste0(
            'curl --data "confirm=1" --output ',
            '~/dva/files/t2d_diagram.zip',
            ' https://www.diagram-consortium.org/check_check.scott_1.php'
        )
    )
    unzip("~/dva/files/t2d_diagram.zip", exdir = "~/dva/files/")
    unlink("~/dva/files/t2d_diagram.zip")
    file.rename(
        "~/dva/files/METAANALYSIS_DIAGRAM_SE1.txt",
        "~/dva/files/t2d_diagram.txt"
    )
}    

## Reference panel for clumpling: 1KG without palindromic or INDELs
if(!dir.exists("~/dva/files/1kg_ref"))
    dir.create("~/dva/files/1kg_ref")

if(!file.exists("~/dva/files/1kg.tgz"))
    download.file("http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz",
                  "~/dva/files/1kg.tgz")

if(!identical(list.files("~/dva/files/1kg_ref"),
              c("EUR.bed", "EUR.bim", "EUR.fam")))
    untar("~/dva/files/1kg.tgz",
          files = c("EUR.bed", "EUR.bim", "EUR.fam"),
          exdir = "~/dva/files/1kg_ref/")

