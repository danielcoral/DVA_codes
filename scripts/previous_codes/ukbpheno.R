# Using ukbtools to extract UKBB data

library(tidyverse)
library(ukbtools)

## Pascal's application

ukbapp1 <- "~/dva/files/ukb_57232"

baskets1 <- gsub("\\.r$", "", list.files(ukbapp1, pattern = "\\.r$"))

ukbphen1 <- ukb_df(baskets1, path = ukbapp1)

## Additional data from Naeimeh's application

ukbapp2 <- "~/dva/files/ukb_18274"

baskets2 <- gsub("\\.r$", "", list.files(ukbapp2, pattern = "\\.r$"))

ukbphen2 <- lapply(baskets2, ukb_df, path = ukbapp2) %>%
    reduce(full_join, by = "eid")

## Mapping eids to Pascal's application
ids <- bind_cols(
    rio::import("~/dva/files/ukb_57232/ukb57232_cal_chr1_v2_s488244.fam",
                format = "\t", select = c(1,5), col.names = c("id1", "sex1")),
    rio::import("~/dva/files/ukb_18274/ukb18274_cal_chr1_v2_s488234.fam",
                format = "\t", select = c(1,5), col.names = c("id2", "sex2"))
) %>%
    ## Ensuring same sex in both datasets
    filter(id1 > 0, id2 > 0, sex1 == sex2)

## Joining datasets
ukbphen <- inner_join(ukbphen1, ids, by = c("eid" = "id1")) %>%
    inner_join(ukbphen2, by = c("id2" = "eid"))

save(ukbphen, file = "~/dva/files/ukbphen.RData")


