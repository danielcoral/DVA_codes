## Recoding variables of interest

## Libraries
library(tidyverse)

## Cleaned UKB dataset
ukb_c <- rio::import("../files/ukb_c.RData")

## Treatments codes
ht_codes <- rio::import("../files/ht_codes.txt")
chol_codes <- rio::import("../files/chol_codes.txt")
t2d_codes <- rio::import("../files/t2d_codes.txt")

## ICD10 Codes
ht_icd10 <- "I10"
t2d_icd10 <- "E11"
t1d_icd10 <- "E10"
chd_icd10 <- paste(paste0("I", 20:25), collapse = "|")

## Columns with information on ICD10 diagnoses
icd_dx <- which(grepl("^diagnoses_icd10", names(ukb_c)))
icd_dates <- which(grepl("^date_.*icd10", names(ukb_c)))

## Columns with information on medication
sr_med_cols <- which(grepl("^medication_for.*_0_[0-9]+$", names(ukb_c)))
med_codes <- which(grepl("treatmentmedication.*_0_[0-9]+$", names(ukb_c)))

## Columns with information on vascular or heart problems
vascheart <- which(grepl("^vascularheart_problems_diagnosed_by_doctor_f6150_0",
                         names(ukb_c)))

## Columns with information on self-reported non-cancer illness codes
sr_ill_codes <- which(grepl("^noncancer.*_0_[0-9]+$", names(ukb_c))) 

## Recoding
ukb_f <- ukb_c %>%
    rename(age0 = age_at_recruitment_f21022_0_0,
           date0 = date_of_attending_assessment_centre_f53_0_0) %>%
    mutate(across(all_of(icd_dates), function(x) age0 - as.numeric((date0 - x) / 365.25),
                  .names = "{.col}_calc_age")) %>%
    mutate(across(all_of(c(icd_dx, icd_dates, sr_med_cols,
                           med_codes, vascheart, sr_ill_codes)),
                  as.character)) %>%
    rowwise() %>%
    mutate(
        sr_bp_med = any(c_across(all_of(sr_med_cols)) == "Blood pressure medication",
                        na.rm = T),
        cod_bp_med = any(c_across(all_of(med_codes)) %in%
                         as.character(ht_codes$coding), na.rm = T),
        icd_ht = any(grepl(ht_icd10, c_across(all_of(icd_dx))) &
                     as.Date(c_across(all_of(icd_dates))) < date0, na.rm = T),
        hbp_doc = any(c_across(all_of(vascheart)) == "High blood pressure", na.rm = T),
        hbp_sr = any(c_across(all_of(sr_ill_codes)) %in% c("1065", "1072"), na.rm = T),
        chol_med1 = any(c_across(all_of(sr_med_cols)) == "Cholesterol lowering medication",
                        na.rm = T),
        chol_med2 = any(c_across(all_of(med_codes)) %in%
                        as.character(chol_codes$coding), na.rm = T),
        icd_t2d = any(grepl(t2d_icd10, c_across(all_of(icd_dx))) &
                      as.Date(c_across(all_of(icd_dates))) < date0 &
                      c_across(ends_with("_calc_age")) > 25),
        t2d_sr = any(c_across(all_of(sr_ill_codes)) == "1223", na.rm = T),
        cur_ins = any(c_across(all_of(sr_med_cols)) == "Insulin", na.rm = T),
        cur_hypog = any(c_across(all_of(med_codes)) %in%
                        as.character(t2d_codes$coding), na.rm = T),
        icd_t1d = any(grepl(t1d_icd10, c_across(all_of(icd_dx))) &
                      as.Date(c_across(all_of(icd_dates))) < date0, na.rm = T),
        t1d_sr = any(c_across(all_of(sr_ill_codes)) == "1222", na.rm = T),
        icd_chd = any(grepl(chd_icd10, c_across(all_of(icd_dx))) &
                      as.Date(c_across(all_of(icd_dates))) < date0, na.rm = T),
        chd_doc = any(c_across(all_of(vascheart)) %in% c("Angina", "Heart attack"),
                      na.rm = T),
        chd_sr = any(c_across(all_of(sr_ill_codes)) %in% c("1074", "1075"), na.rm = T)
    ) %>%
    ungroup() %>%
    transmute(
        eid, age0, date0,
        sex = factor(sex_f31_0_0, ordered = F),
        bmi = body_mass_index_bmi_f21001_0_0,
        obesity = factor(bmi >= 30),
        whr = waist_circumference_f48_0_0 / hip_circumference_f49_0_0,
        armfat = rowMeans(cbind(arm_fat_mass_right_f23120_0_0,
                                arm_fat_mass_left_f23124_0_0)),
        legfat = rowMeans(cbind(leg_fat_mass_right_f23112_0_0,
                                leg_fat_mass_left_f23116_0_0)),
        bp_med = sr_bp_med | cod_bp_med,
        sbp0_0 = coalesce(systolic_blood_pressure_automated_reading_f4080_0_0,
                          systolic_blood_pressure_manual_reading_f93_0_0),
        sbp0_1 = coalesce(systolic_blood_pressure_automated_reading_f4080_0_1,
                          systolic_blood_pressure_manual_reading_f93_0_1),
        sbp = rowMeans(cbind(sbp0_0, sbp0_1)),
        sbp = ifelse(bp_med, sbp + 15, sbp),
        dbp0_0 = coalesce(diastolic_blood_pressure_automated_reading_f4079_0_0,
                          diastolic_blood_pressure_manual_reading_f94_0_0),
        dbp0_1 = coalesce(diastolic_blood_pressure_automated_reading_f4079_0_1,
                          diastolic_blood_pressure_manual_reading_f94_0_1),
        dbp = rowMeans(cbind(dbp0_0, dbp0_1)),
        dbp = ifelse(bp_med, dbp + 10, dbp),
        ## Single visit diagnosis of HBP is also considered only if Grade 3 HT
        hbp = icd_ht | hbp_doc | hbp_sr | sr_bp_med | cod_bp_med | sbp > 180 | dbp > 110,
        hbp = factor(replace_na(hbp, F)),
        hdl = hdl_cholesterol_f30760_0_0,
        tg = triglycerides_f30870_0_0,
        chol_med = chol_med1 | chol_med2,
        t_chol = ifelse(chol_med, cholesterol_f30690_0_0 / 0.8, cholesterol_f30690_0_0),
        ldl = ifelse(chol_med, ldl_direct_f30780_0_0 / 0.7, ldl_direct_f30780_0_0),
        dm_doc = diabetes_diagnosed_by_doctor_f2443_0_0 == "Yes",
        gdm = gestational_diabetes_only_f4041_0_0 == "Yes",
        age_dmdx = coalesce(age_diabetes_diagnosed_f2976_0_0,
                            age_diabetes_diagnosed_f2976_1_0,
                            age_diabetes_diagnosed_f2976_2_0),
        early_dmdx = age_dmdx < 25, 
        dmdx_bf0 = age_dmdx < age0,
        h_glyc = glucose_f30740_0_0 >= 11,
        h_hba1c = glycated_haemoglobin_hba1c_f30750_0_0 >= 48,
        t2d_verbal = (dm_doc | dmdx_bf0) & !(early_dmdx | gdm),
        t2d_meds = cur_ins | cur_hypog,
        t2d_measured = h_glyc | h_hba1c,
        t1d = icd_t1d | t1d_sr,
        t2d = factor(
            replace_na((t2d_verbal | t2d_sr | icd_t2d | t2d_meds | t2d_measured) & !t1d, F)
        ),
        chd = factor(replace_na(icd_chd | chd_doc | chd_sr, F)),
        genetic_principal_components_f22009_0_1,
        genetic_principal_components_f22009_0_2,
        genetic_principal_components_f22009_0_3,
        genetic_principal_components_f22009_0_4,
        genetic_principal_components_f22009_0_5,
        genetic_principal_components_f22009_0_6,
        genetic_principal_components_f22009_0_7,
        genetic_principal_components_f22009_0_8,
        genetic_principal_components_f22009_0_9,
        genetic_principal_components_f22009_0_10
    ) %>%
    select(-c(sbp0_0, sbp0_1, dbp0_0, dbp0_1, bp_med,
              chol_med, dm_doc, gdm, age_dmdx,
              early_dmdx, dmdx_bf0, h_glyc, h_hba1c,
              t2d_verbal, t2d_meds, t2d_measured, t1d))

save(ukb_f, file = "../files/ukb_f.RData")
