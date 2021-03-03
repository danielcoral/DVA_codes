# Using ukbtools to extract UKBB data

library(ukbtools)

ukbphen <- ukb_df("ukb41806", path = "../files/ukb_57232")

save(ukbphen, file = "../files/ukbphen.RData")

df_field <- ukb_df_field("ukb41806", path = "../files/ukb_57232")

save(df_field, file = "../files/df_field.RData")
