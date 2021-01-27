#### UKBB cleaning

library(tidyverse)

### Phenotype data
ukbphen <- rio::import("../files/ukbphen.RData")

## QC of samples
qc_pass <- c(
    ukbphen$sex_f31_0_0 == ukbphen$genetic_sex_f22001_0_0 &
    ukbphen$genetic_ethnic_grouping_f22006_0_0 == "Caucasian" &
    ukbphen$genetic_kinship_to_other_participants_f22021_0_0 == "No kinship found" &
    ukbphen$used_in_genetic_principal_components_f22020_0_0 == "Yes" &
    is.na(ukbphen$outliers_for_heterozygosity_or_missing_rate_f22027_0_0) &
    is.na(ukbphen$sex_chromosome_aneuploidy_f22019_0_0) &
    is.na(ukbphen$country_of_birth_nonuk_origin_f20115_0_0)
)

ukb_c <- ukbphen[which(qc_pass),]

### Function to remove variables without any data (all NA)
remove_all_na <- function(dat){
    col_to_r <- sapply(1:ncol(dat), function(x) all(is.na(dat[,x])))
    if(any(col_to_r)){
        return(dat[,-which(col_to_r)])
    } else {
        return(dat)
    }
}

ukb_c <- remove_all_na(ukb_c)

### Identifying NAs in numeric variables

## NA coded as negative values
ukb_c <- ukb_c %>%
    mutate(
        across(
            ## Only genetic principal components can have negative values
            where(is.numeric) & !starts_with("genetic_principal_components"),
            function(x) ifelse(x < 0, NA, x)
        )
    )

### Outliers using IQR with factor of 2.22 (~3SD)
to_clean <- select(
    ukb_c,
    where(is.numeric) & !(matches("eid|age") | starts_with("genetic_principal_components"))
) %>%
    names() %>%
    str_remove_all("_f[0-9]*_[0-9]_[0-9]$") %>%
    unique()

for(t in to_clean){
    ind <- grep(paste0("^", t, "_f[0-9]*_[0-9]_[0-9]$"), colnames(ukb_c))
    p25 <- quantile(as.matrix(ukb_c[,ind]), 0.25, na.rm = T)
    p75 <- quantile(as.matrix(ukb_c[,ind]), 0.75, na.rm = T)
    iqr <- p75 - p25
    limit <- 2.22 * iqr
    upper <- p75 + limit
    lower <- ifelse(p25 - limit < 0, 0, p25 - limit)
    ukb_c <- ukb_c %>%
        mutate(across(all_of(ind), function(x) ifelse(x < lower | x > upper, NA, x)))
}

## Removing again uninformative variables
ukb_c <- remove_all_na(ukb_c)

save(ukb_c, file = "../files/ukb_c.RData")
