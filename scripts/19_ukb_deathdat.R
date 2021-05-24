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

mor <- rio::import("~/dva/files/death.txt")
morc <- rio::import("~/dva/files/death_cause.txt")
mor_ed <- mor %>%
    ## Only primary cause of death
    filter(ins_index == 0) %>%
    transmute(eid,
              date_of_death = as.Date(format(as.Date(date_of_death,
                                                     format = "%d/%m/%Y"),
                                             "%Y-%m-%d")),
              ## Primary cause of death: Cardiovascular
              cv_m = eid %in% unique(morc$eid[grep("I", morc$cause_icd10)]),
              ## Primary cause of death: Diabetes complications
              t2d_m = eid %in% unique(morc$eid[grep("E11", morc$cause_icd10)]))

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
              t2d_cause = ifelse(death_recorded & death_bf_limit & t2d_m, 1, 0),
              time = round(ifelse(all_cause == 1,
                                  as.numeric((date_of_death - date0)),
                                  as.numeric((as.Date("2021-02-28") - date0))))) %>%
    rename_with(~ gsub("genetic_principal_components_f22009_0_", "gpc", .x)) %>%
    ## GRS
    inner_join(all_scores)

save(ukbphen, file = "~/dva/files/ukbphen.RData")
