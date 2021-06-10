## Install packages needed
install.packages(
    c("tidyverse", "ggrepel", "data.table", "vroom", "rio",
      "patchwork", "ukbtools", "BiocManager",
      "survival", "meta", "devtools", "googledrive",
      "easyPubMed", "pals", "ggforce")
)

BiocManager::install(c("snpStats", "rtracklayer", "biomaRt", "ggtree"))

devtools::install_github(
              c("robingenuer/CoVVSURF",
                "explodecomputer/genetics.binaRies",
                "mrcieu/ieugwasr",
                "teunbrand/ggh4x",
                "ramiromagno/gwasrapidd")
          )

rio::install_formats()
