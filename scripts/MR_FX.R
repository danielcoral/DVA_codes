#### SMR & HEIDI Functions

### SMR p-value given the effects of variants on exposure x and outcome y
pval_smr <- function(b_zx, se_zx, b_zy, se_zy){
    ## Z - scores squared
    z2_zy <- (b_zy / se_zy)^2
    z2_zx <- (b_zx / se_zx)^2
    ## T statistic
    t <- (z2_zy * z2_zx) / (z2_zy + z2_zx)
    ## P value based on Chi squared distribution
    ## 1 degree of freedom
    p <- pchisq(t, df = 1, lower.tail = F)
    return(p)
}

### SMR standard error of beta
std_err_smr <- function(b_zx, se_zx, b_zy, se_zy){
    ## SMR beta
    b_xy <- b_zy / b_zx
    ## Squared coefficients of variation
    varcoef_zx <- (se_zx/b_zx)^2
    varcoef_zy <- (se_zy/b_zy)^2
    ## Sum of coefficients of variation
    sum_varcoef <- varcoef_zx + varcoef_zy
    ## Variance
    var_xy <- (b_xy^2) * sum_varcoef
    se <- sqrt(var_xy)
    return(se)
}

### Covariance
cov_bXY <- function(bxy_i, bxy_j, bzx, bzx_se, bzy, bzy_se, ldrho) {
    zscoreZX <- bzx / bzx_se
    nsnp <- dim(ldrho)[1]
    covbXY <- diag(nsnp)
    if(nsnp>1) {
        zszxij <- zscoreZX%*%t(zscoreZX)
        bxyij <- bxy_i*bxy_j
        sezyij <- bzy_se%*%t(bzy_se)
        bzxij <- bzx%*%t(bzx)
        covbXY <- ldrho*sezyij/bzxij + ldrho*bzxij/zszxij
    }
    return(covbXY)
}

### HEIDI p-value
heidi_pvalue <- function(bzx, bzx_se, bzy, bzy_se, ldrho, topsnp_index) {
    ## 1 the top SNP
    bxy_top <- bzy[topsnp_index]/bzx[topsnp_index]
    bxy_top_se <- std_err_smr(bzx[topsnp_index], bzx_se[topsnp_index],
                              bzy[topsnp_index], bzy_se[topsnp_index])
    ## 2 d = bxy_i - bxy_top
    bxy_snp <- bzy / bzx
    d <- bxy_snp - bxy_top
    cov_bxy <- cov_bXY(bxy_snp, bxy_snp, bzx, bzx_se, bzy, bzy_se, ldrho)
    var_d <- cov_bxy - cov_bxy[topsnp_index,] + bxy_top_se^2
    var_d <- t(t(var_d) - cov_bxy[topsnp_index,])
    ## 3 p-value for top SNP
    var_d[topsnp_index, topsnp_index] <- 1
    chival <- d^2/diag(var_d)
    ## 4. Saddlepoint method
    m <- length(bzx) - 1
    d <- d[-topsnp_index]
    var_d <- var_d[-topsnp_index, -topsnp_index]
    chival <- chival[-topsnp_index]
    ## Eigen decomposition
    corr_d <- diag(m)
    for( i in 1 : (m-1) ) {
        for( j in (i+1) : m ) {
            corr_d[i,j] = corr_d[j,i] = var_d[i,j] / sqrt(var_d[i,i]*var_d[j,j])
        }         
    }
    ## estimate the p-value
    lambda <- eigen(corr_d, symmetric=TRUE, only.values=TRUE)$values
    t <- length(lambda)
    heidi_p <- survey::pchisqsum(sum(chival), df=rep(1,t), a=lambda,
                                 method="saddlepoint", lower.tail = F)
    return(heidi_p)   
}

library(MendelianRandomization)

ivwfx <- function(b_exp, b_out, se_exp, se_out){
    mod <- lm(b_out ~ -1 + b_exp, weights = 1/se_out^2)
    ivw.res <- summary(mod)
    data.frame(method = "mr_ivw",
               b = ivw.res$coef["b_exp","Estimate"],
               se = ivw.res$coef["b_exp","Std. Error"]/min(1,ivw.res$sigma),
               Q_df = length(b_exp) - 1,
               nsnp = length(b_exp)) %>%
        mutate(pval = 2 * pnorm(abs(b/se), lower.tail=FALSE),
               Q = ivw.res$sigma^2 * Q_df,
               Q_pval = pchisq(Q, Q_df, lower.tail=FALSE),
               ci_lo = b - qnorm(1-0.05/2)*se,
               ci_up = b + qnorm(1-0.05/2)*se)
}

eggerfx <- function(b_exp, b_out, se_exp, se_out){
    mod <- lm(b_out ~ b_exp, weights=1/se_out^2)
    smod <- summary(mod)
    data.frame(method = c("intercept", "egger"),
               b = coefficients(smod)[,1],
               se = coefficients(smod)[,2] / min(1,smod$sigma),
               Q_df = length(b_exp) - 2,
               nsnp = length(b_exp)) %>%
        mutate(pval = 2 * pnorm(abs(b / se), lower.tail = FALSE),
               Q = smod$sigma^2 * Q_df,
               Q_pval = pchisq(Q, Q_df, lower.tail=FALSE),
               ci_lo = b - qnorm(1-0.05/2)*se,
               ci_up = b + qnorm(1-0.05/2)*se)
}

wmedianfx <- function(b_exp, b_out, se_exp, se_out){
    b_iv <- b_out / b_exp
    VBj <- ((se_out)^2)/(b_exp)^2 + (b_out^2)*((se_exp^2))/(b_exp)^4
    data.frame(method = "wmedian",
               b = TwoSampleMR::weighted_median(b_iv, 1 / VBj),
               se = TwoSampleMR::weighted_median_bootstrap(b_exp, b_out, se_exp, se_out,
                                                           1 / VBj, 1000),
               nsnp = length(b_exp)) %>%
        mutate(pval = 2 * pnorm(abs(b/se), lower.tail=FALSE),
               ci_lo = b - qnorm(1-0.05/2)*se,
               ci_up = b + qnorm(1-0.05/2)*se)
}

conmixfx <- function(b_exp, b_out, se_exp, se_out){
    input <- mr_input(bx = b_exp, bxse = se_exp,
                      by = b_out, byse = se_out)
    psirange <- seq(1, 2, .1)
    map_dfr(psirange,
            ~ {
                psival <- .x * sd(b_out/b_exp)
                res <- mr_conmix(input, psi = psival)
                data.frame(b = res@Estimate,
                           ci_lo = res@CILower, ci_up = res@CIUpper,
                           pval = res@Pvalue, nsnp = res@SNPs, psi = .x,
                           psitrue = res@Psi, method = "conmix")
            }) %>%
        mutate(anymulti = na_if(any(duplicated(psi)), FALSE))
}

Fstat <- function(b_exp, se_exp) mean((b_exp^2) / (se_exp^2))

isq_gx <- function(b_exp, se_exp){
    k <- length(b_exp)
    w <- 1/se_exp^2
    sum.w <- sum(w)
    mu.hat <- sum(b_exp*w)/sum.w
    Q <- sum(w*(b_exp-mu.hat)^2)
    Isq <- (Q - (k-1))/Q
    Isq <- max(0,Isq)
    return(Isq)
}

mr_analysis <- function(b_exp, b_out, se_exp, se_out){
    bind_rows(ivwfx(b_exp, b_out, se_exp, se_out),
              eggerfx(b_exp, b_out, se_exp, se_out),
              wmedianfx(b_exp, b_out, se_exp, se_out),
              conmixfx(b_exp, b_out, se_exp, se_out)) %>%
        mutate(fstat = Fstat(b_exp, se_exp),
               isq = isq_gx(b_exp, se_exp))
}

## Multivariable MR

mvmr.ivwfx <- function(Bx, By, Bxse, Byse) {
    summary <- summary(lm(By ~ Bx - 1, weights = Byse^(-2)))
    data.frame(method = "mr_ivw",
               b = summary$coef[,1],
               se = summary$coef[,2]/min(summary$sigma,1),
               nsnp = dim(Bx)[1]) %>%
        mutate(pval = 2*pnorm(-abs(b/se)),
               Q_df = dim(Bx)[1] - dim(Bx)[2] - 1,
               Q = Q_df * (summary$sigma^2),
               Q_pval = pchisq(Q, df = Q_df, lower.tail = FALSE),
               ci_lo = b - qnorm(1-0.05/2) * se,
               ci_up = b + qnorm(1-0.05/2)* se)
}

mvmr.eggerfx <- function(Bx, By, Bxse, Byse) {
    summary <- summary(lm(By ~ Bx, weights = Byse^(-2)))
    data.frame(method = c("intercept",
                          rep("egger", dim(Bx)[2])),
               b = summary$coef[,1],
               se = summary$coef[,2]/min(summary$sigma,1),
               nsnp = dim(Bx)[1]) %>%
        mutate(pval = 2*pnorm(-abs(b/se)),
               Q_df = dim(Bx)[1] - dim(Bx)[2] - 1,
               Q = Q_df * (summary$sigma^2),
               Q_pval = pchisq(Q, df = Q_df, lower.tail = FALSE),
               ci_lo = b - qnorm(1-0.05/2) * se,
               ci_up = b + qnorm(1-0.05/2)* se)
}

mvmr.wmedianfx <- function(Bx, By, Bxse, Byse){
    i <- mr_mvinput(bx = Bx, bxse = Bxse, by = By, byse = Byse)
    resmed <- mr_mvmedian(i)
    data.frame(method = "wmedian",
               b = resmed@Estimate, se = resmed@StdError,
               nsnp = resmed@SNPs,
               pval = resmed@Pvalue,
               ci_lo = resmed@CILower, ci_up = resmed@CIUpper)
}

mvmr_analysis <- function(Bx, By, Bxse, Byse, exposures = NULL){
    ivw <- mvmr.ivwfx(Bx, By, Bxse, Byse)
    ivw$exposures <- exposures
    egger <- mvmr.eggerfx(Bx, By, Bxse, Byse)
    egger$exposures <- c("intercept", exposures)
    wmed <- mvmr.wmedianfx(Bx, By, Bxse, Byse)
    wmed$exposures <- exposures
    bind_rows(ivw, egger, wmed)
}
