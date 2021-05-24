## Regression of BMI score on T2D by trait quintiles

library(tidyverse)

multiprs_res <- lapply(
    list.files("../files/multiprs/", full.names = TRUE),
    function(x){
        traitname <- gsub("../files/multiprs//|\\.tsv", "", x)
        rio::import(x) %>%
            `names<-`(c("eid", traitname))
    }
) %>%
    reduce(inner_join) %>%
    mutate(across(-eid, scale))

multiprs_res <- multiprs_res %>% select(-whr_prs)

head(multiprs_res)

reg_data <- rio::import("../files/reg_data.tsv")

dat <- inner_join(reg_data, multiprs_res, by = "eid")

head(dat, n = 2)

res <- glm(t2d_dx0 ~ sbp_prs + bmi_prs + sbp_prs * bmi_prs +
               age0 + sex +
               genpcs_1 + genpcs_2 + genpcs_3 + genpcs_4 + genpcs_5 +
               genpcs_6 + genpcs_7 + genpcs_8 + genpcs_9 + genpcs_10,
           data = dat, family = "binomial")

res <- dat %>%
    mutate("sbpquant" = ntile(sbp_prs, 2)) %>%
    group_by("sbpquant") %>%
    group_map(
        ~ glm(t2d_dx0 ~ bmi_prs + age0 + sex +
                  genpcs_1 + genpcs_2 + genpcs_3 + genpcs_4 + genpcs_5 +
                  genpcs_6 + genpcs_7 + genpcs_8 + genpcs_9 + genpcs_10,
              data = .x, family = "binomial")
    )

summary(res)
