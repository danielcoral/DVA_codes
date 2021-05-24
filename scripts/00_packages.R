## Install packages needed
install.packages(
    c("tidyverse", "ggrepel", "data.table", "vroom", "rio", "VennDiagram",
      "patchwork", "survey", "dendextend", "ukbtools", "BiocManager",
      "survival", "ggfortify", "meta", "devtools", "pROC", "googledrive",
      "easyPubMed", "pals", "ggforce")
)

BiocManager::install(c("snpStats", "rtracklayer", "biomaRt", "ggtree"))

devtools::install_github(
              c("phenoscanner/phenoscanner",
                "robingenuer/CoVVSURF",
                "ricardo-bion/ggradar",
                "explodecomputer/genetics.binaRies",
                "mrcieu/ieugwasr",
                "mrcieu/gwasvcf",
                "mrcieu/gwasglue",
                "anastasia-lucas/hudson",
                "teunbrand/ggh4x",
                "ramiromagno/gwasrapidd")
          )

rio::install_formats()
