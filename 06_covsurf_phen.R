#### CoVVSURF - Manually filtered phenome matrix

### The matrix
mat_phen <- rio::import("../files/mat_phen.txt")

### Obtaining the hierarchical clustering tree
cov_phen <- CoVVSURF::CoV(mat_phen[,-c(1,2)])

save(cov_phen, file = "../files/cov_phen.RData")

### Initial search for best range (for optimization)
kval_init <- c(2, seq(10, ncol(mat_phen) - 2, 10))

res_phen_init <- CoVVSURF::covsurf(X = mat_phen[,-c(1,2)],
                                   y = as.factor(mat_phen$disc),
                                   kval = kval_init,
                                   tree = cov_phen)

save(res_phen_init, file = "../files/res_phen_init.RData")

kopt_index <- which(kval_init == res_phen_init$kopt)

bestrange <- kval_init[c(max(1, kopt_index - 2),
                         min(length(kval_init), kopt_index + 2))]

### Running COVVSURF across the best range found
res_phen <- CoVVSURF::covsurf(X = mat_phen[,-c(1,2)],
                              y = as.factor(mat_phen$disc),
                              kval = bestrange[1]:bestrange[2],
                              tree = cov_phen)

save(res_phen, file = "../files/res_phen.RData")
