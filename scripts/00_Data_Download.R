#### Data download

## BMI from the GIANT Consortium + UK Biobank
download.file(
  "https://portals.broadinstitute.org/collaboration/giant/images/c/c8/Meta-analysis_Locke_et_al%2BUKBiobank_2018_UPDATED.txt.gz",
  "/ludc/Home/daniel_c/projects/DVA/Data/GWAS_sumstats/bmi.txt.gz"
)
R.utils::gunzip("/ludc/Home/daniel_c/projects/DVA/Data/GWAS_sumstats/bmi.txt.gz")

## T2D from DIAGRAM Consortium + UK Biobank
system(
  paste0(
    'curl --data "confirm=1" --output ',
    '/ludc/Home/daniel_c/projects/DVA/Data/GWAS_sumstats/t2d.zip',
    ' https://www.diagram-consortium.org/check_check.Mahajan2018b.T2D.3.php'
  )
)
unzip(
  "/ludc/Home/daniel_c/projects/DVA/Data/GWAS_sumstats/t2d.zip", 
  exdir = "/ludc/Home/daniel_c/projects/DVA/Data/GWAS_sumstats"
)
unlink("/ludc/Home/daniel_c/projects/DVA/Data/GWAS_sumstats/t2d.zip")
file.rename(
  "/ludc/Home/daniel_c/projects/DVA/Data/GWAS_sumstats/Mahajan.NatGenet2018b.T2D.European.txt",
  "/ludc/Home/daniel_c/projects/DVA/Data/GWAS_sumstats/t2d.txt"
)

## T2D from DIAGRAM only
system(
  paste0(
    'curl --data "confirm=1" --output ',
    '/ludc/Home/daniel_c/projects/DVA/Data/GWAS_sumstats/t2d_diagram.zip',
    ' https://www.diagram-consortium.org/check_check.scott_1.php'
  )
)
unzip("/ludc/Home/daniel_c/projects/DVA/Data/GWAS_sumstats/t2d_diagram.zip", exdir = "~/dva/files/")
unlink("/ludc/Home/daniel_c/projects/DVA/Data/GWAS_sumstats/t2d_diagram.zip")
file.rename(
  "/ludc/Home/daniel_c/projects/DVA/Data/GWAS_sumstats/METAANALYSIS_DIAGRAM_SE1.txt",
  "/ludc/Home/daniel_c/projects/DVA/Data/GWAS_sumstats/t2d_diagram.txt"
)    

## Reference panel for clumpling: 1KG without palindromic or INDELs
dir.create("/ludc/Home/daniel_c/projects/DVA/Data/ReferenceData/1kg_ref")
download.file("http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz", "projects/DVA/Data/ReferenceData/1kg_ref/1kg.tgz")
untar(
  "projects/DVA/Data/ReferenceData/1kg_ref/1kg.tgz", 
  files = c("EUR.bed", "EUR.bim", "EUR.fam"),
  exdir = "projects/DVA/Data/ReferenceData/1kg_ref/"
)

## Gencode
download.file(
  "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz",
  "~/projects/DVA/Data/ReferenceData/gencode.v19.annotation.gtf.gz"
)
