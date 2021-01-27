#### CoVVSURF - Manually filtered phenome matrix

### The matrix
mat_phen <- rio::import("../files/mat_phen.txt")

### Clustering
cov_phen <- CoVVSURF::CoV(mat_phen[,-c(1,2)])

save(cov_phen, file = "../files/cov_phen.RData")

### Discordance classification
res_phen <- CoVVSURF::covsurf(mat_phen[,-c(1,2)],
                            y = as.factor(mat_phen$disc),
                            tree = cov_phen)

save(res_phen, file = "../files/res_phen.RData")
