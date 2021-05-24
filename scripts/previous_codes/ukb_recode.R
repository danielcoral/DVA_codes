## Recoding variables of interest

library(tidyverse)

ukb_c <- rio::import("~/dva/files/ukb_c.RData")
t2d_diagnosis <- rio::import("~/dva/files/t2d_diagnosis.tsv")

## Treatments codes
ht_codes <- rio::import("~/dva/files/ht_codes.txt") %>% mutate(coding = as.character(coding))
chol_codes <- rio::import("~/dva/files/chol_codes.txt") %>% mutate(coding = as.character(coding))

## Columns with information on medication
ts.meds.colnam <- which(grepl("^medication_for.*_0_[0-9]+$", names(ukb_c)))
vi.meds.colnam <- which(grepl("treatmentmedication.*_0_[0-9]+$", names(ukb_c)))

ts.meds.cols <- ukb_c %>%
    transmute(eid, across(all_of(ts.meds.colnam), as.character)) %>%
    pivot_longer(-eid, names_to = NULL, values_to = "meds")

vi.meds.cols <- ukb_c %>%
    transmute(eid, across(all_of(vi.meds.colnam), as.character)) %>%
    pivot_longer(-eid, names_to = NULL, values_to = "meds")

medstab <- list(
    ## Blood pressure treatment - Touchscreen
    ts.meds.cols %>%
    group_by(eid) %>%
    summarise(ts.bpmed = any(meds == "Blood pressure medication", na.rm = TRUE),
              ts.cholmed = any(meds == "Cholesterol lowering medication", na.rm = TRUE)),
    vi.meds.cols %>%
    group_by(eid) %>%
    summarise(vi.bpmed = any(meds %in% ht_codes$coding, na.rm = TRUE),
              vi.cholmed = any(meds %in% chol_codes$coding, na.rm = TRUE))
) %>%
    modify_at(2, ~ select(.x, -eid)) %>%
    bind_cols()

## Recoding
ukb_f <- ukb_c %>%
    inner_join(medstab, by = "eid") %>%
    inner_join(t2d_diagnosis, by = "eid") %>%
    rename_with(~ gsub(".*f22009_0", "genpcs", .x),
                matches("f22009")) %>%
    transmute(eid,
              sex = factor(sex_f31_0_0, ordered = F),
              age0 = age_at_recruitment_f21022_0_0,
              date0 = date_of_attending_assessment_centre_f53_0_0,
              across(all_of(paste0("genpcs_", 1:10))),
              t2d_dx0 = t2d_dx,
              bmi0 = body_mass_index_bmi_f21001_0_0,
              whr0 = waist_circumference_f48_0_0 / hip_circumference_f49_0_0,
              armfat0 = rowMeans(cbind(arm_fat_mass_right_f23120_0_0,
                                      arm_fat_mass_left_f23124_0_0)),
              legfat0 = rowMeans(cbind(leg_fat_mass_right_f23112_0_0,
                                      leg_fat_mass_left_f23116_0_0)),
              sbp0 = rowMeans(
                  cbind(
                      coalesce(systolic_blood_pressure_automated_reading_f4080_0_0,
                               systolic_blood_pressure_manual_reading_f93_0_0),
                      coalesce(systolic_blood_pressure_automated_reading_f4080_0_1,
                               systolic_blood_pressure_manual_reading_f93_0_1)
                  )
              ),
              sbp0 = ifelse(ts.bpmed | vi.bpmed, sbp0 + 15, sbp0),
              dbp0 = rowMeans(
                  cbind(
                      coalesce(diastolic_blood_pressure_automated_reading_f4079_0_0,
                               diastolic_blood_pressure_manual_reading_f94_0_0),
                      coalesce(diastolic_blood_pressure_automated_reading_f4079_0_1,
                               diastolic_blood_pressure_manual_reading_f94_0_1)
                  )
              ),
              dbp0 = ifelse(ts.bpmed | vi.bpmed, dbp0 + 10, dbp0),
              hdl0 = hdl_cholesterol_f30760_0_0,
              tg0 = triglycerides_f30870_0_0,
              ldl0 = ldl_direct_f30780_0_0,
              ldl0 = ifelse(ts.cholmed | vi.cholmed, ldl0 / 0.7, ldl0))

save(ukb_f, file = "~/dva/files/ukb_f.RData")
