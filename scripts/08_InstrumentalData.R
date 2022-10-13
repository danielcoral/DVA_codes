# Data for MR analysis

# Relevant exposures identified in the previous analyses

exposures <- tibble::tribble(
    ~trait_short, ~id, 
    "Urate", "ukb-d-30880_irnt",
    "WHR", "ieu-a-73",
    "#Leuco", "ukb-d-30000_irnt",
    "SBP", "ieu-b-38",
    "DBP", "ieu-b-39",
    "HDL", "ukb-d-30760_irnt",
    "LDL", "ukb-d-30780_irnt",
    "TG", "ukb-d-30870_irnt",
    "HeelBMD", "ukb-b-20124",
    "Phosphate", "ukb-d-30810_irnt",
    "Platelet volume", "ukb-d-30100_irnt",
    "Urea", "ukb-d-30670_irnt",
    "MCV", "ukb-d-30040_irnt",
    "MCH", "ukb-d-30050_irnt",
    "EDW", "ukb-d-30070_irnt",
    "MCHC", "ukb-d-30060_irnt",
    "SHBG", "ukb-d-30830_irnt",
    "GGT", "ukb-d-30730_irnt",
    "ALT", "ukb-d-30620_irnt",
    "ApoA1", "ieu-b-107",
    "PTH", "prot-a-2431",
    "CRP", "ukb-d-30710_irnt",
)

sigins <- ieugwasr::tophits(exposures$id, clump = 0)

readr::write_tsv(sigins, "../data/sigins.tsv")