#### CoVVSURF - Phenome matrix

mat_phen <- rio::import("~/dva/files/mat_phen.tsv")

cov_phen <- CoVVSURF::CoV(mat_phen[-c(1,2)])

save(cov_phen, file = "~/dva/files/cov_phen.RData")

res_phen <- CoVVSURF::covsurf(mat_phen[,-c(1:3)], factor(mat_phen$disc),
                              tree = cov_phen, ncores = parallel::detectCores())

save(res_phen, file = "~/dva/files/res_phen.RData")
