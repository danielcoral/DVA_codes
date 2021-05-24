## Type 2 Diabetes diagnosis in UK Biobank

library(tidyverse)

ukb_c <- rio::import("~/dva/files/ukb_c.RData")

## Diabetes diagnosis can be derived from:
t2d.f <- tribble(
    ~f.nam, ~f.id,
    ## A. The verbal interview:
    "vi.dx", "f20002", ## 1. Reported diagnosis (in codes)
    "vi.dx.age", "f20009", ## 2. Age when diagnosed
    "vi.meds", "f20003", ## 3. Reported medications
    ## B. The touchscreen questionnaire
    "ts.dmdx", "f2443", ## 4. Reported diabetes
    "ts.dmdx.age", "f2976", ## 5. Age of diagnosis of diabetes
    "ts.earlyins", "f2986", ## 6. Started insulin within 1 year diagnosis 
    "ts.gdm.only", "f4041", ## 7. Gestational diabetes only
    "ts.ins.m", "f6177", ## 8. Taking insulin - Men
    "ts.ins.f", "f6153", ## 9. Taking insulin - Women
    ## C. ICD codes
    "icd10", "f41270", ## 10. Code
    "icd10.date", "f41280", ## 11. Date of ICD10 diagnosis
    ## D. Biochemistry
    "glucose", "f30740", ## 12. Glucose
    "fasting", "f74", ## 13. Fasting time
    "hba1c", "f30750" ## 14. HbA1c
) %>%
    ## We will only use first instance
    ## Note: Instances does not apply to ICD10 code data
    mutate(f.id = ifelse(grepl("^icd10", f.nam), f.id,
                         paste0(f.id, "_0_")))

## Type 1 diabetes can be derived from:
## Correspoding code in verbal interview
## Corresponding code in ICD10 codes

## Codes for diabetes medications
diabmeds <- rio::import("~/dva/files/t2d_codes.txt") %>% mutate(coding = as.character(coding))

## Here we used the algorithm proposed by Eastwood et al. 2016

## Verbal interview - Diagnoses
t2dcodes <- c("1220", "1223", "1276", "1486", "1607")

vidx <- t2d.f %>%
    filter(f.nam %in% c("vi.dx", "vi.dx.age")) %>%
    rowwise() %>%
    group_map(~ ukb_c %>%
                  select(eid, matches(.x$f.id)) %>%
                  pivot_longer(-eid, names_to = NULL,
                               values_to = .x$f.nam)) %>%
    modify_at(2, ~ select(.x, -eid)) %>%
    bind_cols() %>%
    group_by(eid) %>%
    summarise(t2d_vi = any(vi.dx %in% t2dcodes, na.rm = TRUE),
              earlyt2d_vi = any(vi.dx %in% t2dcodes & vi.dx.age < 35, na.rm = TRUE),
              t2d_vi_noage = any(vi.dx %in% t2dcodes & is.na(vi.dx.age)),
              t1d = any(vi.dx == "1222", na.rm = TRUE)) %>%
    transmute(eid,
              t2d_vi = case_when(t2d_vi & (earlyt2d_vi | t2d_vi_noage | t1d) ~ NA,
                                 t2d_vi ~ TRUE,
                                 TRUE ~ FALSE))

## Verbal interview - Medications
vimeds <- ukb_c %>%
    select(eid, matches(with(t2d.f, f.id[f.nam == "vi.meds"]))) %>%
    pivot_longer(-eid, names_to = NULL, values_to = "medcodes") %>%
    group_by(eid) %>%
    summarise(t2d_vimeds = any(medcodes %in% diabmeds$coding))

## Touchscreen questionnaire
tsdx.cols <- t2d.f %>%
    filter(grepl("ts\\.[^ins]", f.nam))
tsmed.cols <- t2d.f %>%
    filter(grepl("ts\\.ins", f.nam))

tsdx <- ukb_c %>%
    select(eid, sex_f31_0_0,
           matches(tsdx.cols$f.id)) %>%
    `names<-`(c("eid", "sex", tsdx.cols$f.nam)) %>%
    transmute(eid,
              ## Note: No need to exclude absence of age at diagnosis
              t2d_ts = case_when(
                  ts.dmdx == "Yes" & ts.dmdx.age < 35 | ts.earlyins == "Yes" | ts.gdm.only == "Yes" ~ NA,
                  ts.dmdx == "Yes" & ts.dmdx.age > 35 & ts.earlyins == "No" & ts.gdm.only == "No" ~ TRUE,
                  TRUE ~ FALSE
              ))

tsmeds <- ukb_c %>%
    transmute(eid,
              across(matches(tsmed.cols$f.id), as.character)) %>%
    pivot_longer(-eid, names_to= NULL, values_to = "med") %>%
    group_by(eid) %>%
    summarise(insulin = any(med == "Insulin", na.rm = T))

## ICD-10 Codes
icddx <- t2d.f %>%
    filter(grepl("icd10", f.nam)) %>%
    rowwise() %>%
    group_map(~ ukb_c %>%
                  transmute(eid,
                            age0 = age_at_recruitment_f21022_0_0,
                            date0 = date_of_attending_assessment_centre_f53_0_0,
                            across(matches(.x$f.id))) %>%
                  pivot_longer(matches(.x$f.id),
                               names_to = NULL,
                               values_to = .x$f.nam)) %>%
    modify_at(2, ~ select(.x, -c(eid, age0, date0))) %>%
    bind_cols() %>%
    mutate(icd10.age = age0 - as.numeric((date0 - icd10.date) / 365.25),
           icd10.bf.0 = icd10.age < age0) %>%
    group_by(eid) %>%
    ## Note: All codes have a date of diagnosis
    summarise(icd10.t2d = any(grepl("E11", icd10) &
                              icd10.age > 35 &
                              icd10.bf.0, na.rm = T),
              icd10.t1d = any(grepl("E10", icd10) &
                              icd10.bf.0, na.rm = T)) %>%
    transmute(eid, t2d_icd = ifelse(icd10.t1d, NA, icd10.t2d))

## Biochemistry - Diagnostic criteria by ADA
bioch.cols <- filter(t2d.f,
                     f.nam %in% c("glucose", "fasting", "hba1c"))

biochdx <- ukb_c %>%
    select(eid, matches(bioch.cols$f.id)) %>%
    `names<-`(c("eid", bioch.cols$f.nam)) %>%
    mutate(fastgluc = fasting > 8 & glucose > 7,
           nonfastgluc = fasting < 8 & glucose > 11,
           high_hba1c = hba1c > 48) %>%
    transmute(eid, t2d_bioch = replace_na(fastgluc | nonfastgluc | high_hba1c, FALSE))

## Joining all sources of diagnosis
t2d_diagnosis <- list(vidx, vimeds, tsdx, tsmeds, icddx, biochdx) %>%
    reduce(inner_join) %>%
    pivot_longer(-eid, names_to = NULL, values_to = "dx") %>%
    group_by(eid) %>%
    summarise(t2d_dx = ifelse(anyNA(dx), NA, any(dx)))

## Saving
rio::export(t2d_diagnosis, "~/dva/files/t2d_diagnosis.tsv")
