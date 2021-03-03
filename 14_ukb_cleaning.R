#### UKBB cleaning

library(tidyverse)

### Phenotype data
ukbphen <- rio::import("../files/ukbphen.RData")

ukb_c <- ukbphen %>%
    ## Simple QC filtering - Used in genetic PC calculation (which was QC'd)
    filter(sex_f31_0_0 == genetic_sex_f22001_0_0,
           genetic_ethnic_grouping_f22006_0_0 == "Caucasian",
           used_in_genetic_principal_components_f22020_0_0 == "Yes",
           is.na(sex_chromosome_aneuploidy_f22019_0_0)) %>%
    ## NA coded as negative values
    ## Only genetic principal components can have negative values
    mutate(across(where(is.numeric) & !starts_with("genetic_principal_components"),
                  function(x) ifelse(x < 0, NA, x)))

## Outliers using IQR with factor of 2.22 (~3SD)
to_clean <- select(
    ukb_c,
    where(is.numeric) & !(matches("eid|^age") | starts_with("genetic_principal_components"))
) %>%
    names() %>%
    str_remove_all("_f[0-9]*_[0-9]_[0-9]$") %>%
    unique()

for(t in to_clean){
    ind <- grep(paste0("^", t, "_f[0-9]*_0_[0-9]$"), colnames(ukb_c))
    p25 <- quantile(as.matrix(ukb_c[,ind]), 0.25, na.rm = T)
    p75 <- quantile(as.matrix(ukb_c[,ind]), 0.75, na.rm = T)
    iqr <- p75 - p25
    limit <- 2.22 * iqr
    upper <- p75 + limit
    lower <- ifelse(p25 - limit < 0, 0, p25 - limit)
    ukb_c <- ukb_c %>%
        mutate(across(all_of(ind), function(x) ifelse(x < lower | x > upper, NA, x)))
}

ukb_c <- ukb_c %>%
    ## Only informative variables
    select(eid, where(function(x) !all(is.na(x))))

save(ukb_c, file = "../files/ukb_c.RData")
