# Using ukbtools to extract UKBB data

library(tidyverse)
library(ukbtools)

ukbapp <- "~/dva/files/ukb_57232"

baskets_avail <- gsub("\\.r$", "", list.files(ukbapp, pattern = "\\.r$"))

ukb_covar <- ukb_df(baskets_avail, path = ukbapp)

## PRS
all_scores <- rio::import("~/dva/files/all_scores.tsv") ## %>% mutate(across(-eid, scale))

## UKB death data

## According to UK Biobank guidelines:

## Step 1. Adding field 40023 to an application basket vi the Access Management System (AMS)

##Step 2. Login to the AMS and go to:
##     	      > Projects tab
##	      	> Data tab
##		  > Data refresh or download
##		    > Data portal tab
##		      > Connect
##		      	> Table download
##			  > Once here, type table names: DEATH and DEATH_CAUSE
##			    > Download by clicking on the link

## Columns of DEATH
## eid	         = Individual ID
## ins_index     = Place of death record (individuals can have more than one)
## date_of_death = Self-explanatory
## dsource       = England/Wales or Scotland
## source 	 = Numerical codes for the instances

## Columns of DEATH_CAUSE
## eid	         = Individual ID
## ins_index     = Place of death record (individuals can have more than one)
## arr_index     = 0 for primary cause, > 0 for secondary causes 
## level	 = 1 for primary cause, 2 for contributory causes
## cause_icd10   = Self-explanatory

mordate <- rio::import("~/dva/files/death.txt")
morcause <- rio::import("~/dva/files/death_cause.txt")

## Only first instances and primary cause
mordate <- filter(mordate, ins_index == 0)
morcause <- filter(morcause, ins_index == 0, level == 1)

## Joining
mor <- inner_join(mordate, morcause) %>%
    select(eid, date_of_death, cause_icd10)

## CV death
## Diagnosis from ICD-10 chapter IX (Circulatory system)
cvcauses <- paste0("I", c(1,2,34,35:37,42,44:50,6,70,71,72,739,74,75,77)) %>%
    ## Additional codes extracted from Ruigómez et al. 2021
    c("G45", "G46", "R092", "R960", "R961", "R98", "F01", "R570") %>%
    paste0("^", .) %>%
    paste(collapse = "|")

## Table of codes retrieved for supplementary data
icdcodes <- icd10cm2019

icd_cv <- icdcodes %>%
    filter(grepl(cvcauses, code)) %>%
    select(code, short_desc) %>%
    filter(str_count(code) <= 4)

icd_cv %>%
    rio::export("~/dva/files/icd_cv.tsv")

## Coding date and causes of death
mor_ed <- mor %>%
    transmute(eid,
              date_of_death = as.Date(format(as.Date(date_of_death,
                                                     format = "%d/%m/%Y"),
                                             "%Y-%m-%d")),
              cv_m = grepl(cvcauses, cause_icd10))

## Merging covariate and death data
ukbphen <- ukb_covar %>%
    ## Simple QC filtering - Used in genetic PC calculation (which was QC'd)
    filter(sex_f31_0_0 == genetic_sex_f22001_0_0,
           genetic_ethnic_grouping_f22006_0_0 == "Caucasian",
           used_in_genetic_principal_components_f22020_0_0 == "Yes",
           is.na(sex_chromosome_aneuploidy_f22019_0_0)) %>%
    left_join(mor_ed) %>%
    transmute(eid,
              age0 = age_at_recruitment_f21022_0_0,
              sex = factor(sex_f31_0_0, ordered = FALSE),
              date0 = date_of_attending_assessment_centre_f53_0_0,
              across(starts_with("genetic_principal_components") &
                     matches("_[1-9]$|_10$")),
              death_bf_limit = date_of_death < "2021-02-28", ## Censor date!
              death_recorded = !is.na(date_of_death),
              all_cause = ifelse(death_recorded & death_bf_limit, 1, 0),
              cv_cause = ifelse(all_cause == 1 & cv_m, 1, 0),
              time = round(ifelse(all_cause == 1,
                                  as.numeric((date_of_death - date0)),
                                  as.numeric((as.Date("2021-02-28") - date0))))) %>%
    rename_with(~ gsub("genetic_principal_components_f22009_0_", "gpc", .x)) %>%
    ## GRS
    inner_join(all_scores)

save(ukbphen, file = "~/dva/files/ukbphen.RData")
