## COVVSURF results - Phenome matrix

library(tidyverse)
library(ggtree)
library(patchwork)
library(ggh4x)
library(pROC)
library(NbClust)
library(fpc)


mix <- rio::import("~/dva/files/mix.tsv")
cov_phen <- rio::import("~/dva/files/cov_phen.RData")
res_phen <- rio::import("~/dva/files/res_phen.RData")
mat_phen <- rio::import("~/dva/files/mat_phen.tsv")
premat_phen <- rio::import("~/dva/files/premat_phen.tsv")

##----------------------------
## Optimal number of clusters
##----------------------------

## Labelling threshold
threslabpos <- 175

kopt_phen <- with(res_phen, cbind(kval, oob)) %>%
    data.frame %>%
    `names<-`(c("k", "error", "error_sd")) %>%
    mutate(minerr = error - error_sd, maxerr = error + error_sd) %>%
    {
        kopt_final <- .
        min_err <- min(kopt_final$error)
        min_err.sd <- kopt_final$error_sd[which.min(kopt_final$error)]
        kopt_thresh <- min_err + min_err.sd
        kopt_final_plot <- kopt_final %>%
            ggplot(aes(k, error)) +
            geom_ribbon(aes(ymin = minerr, ymax = maxerr), fill = "grey") +
            geom_line() +
            geom_hline(yintercept = min_err, colour = "green", lty = "dashed") +
            annotate("text", x = threslabpos, y = min_err, label = "Lowest error",
                     vjust = 1.5, hjust = .25) +
            geom_hline(yintercept = kopt_thresh, colour = "maroon", lty = "dashed") +
            annotate("text", x = threslabpos, y = kopt_thresh, label = "Threshold",
                     vjust = -0.5, hjust = .25) +
            geom_segment(aes(x = threslabpos, y = min_err,
                             xend = threslabpos, yend = kopt_thresh),
                         arrow = arrow(ends = "both", length = unit(.2, "cm"))) +
            annotate("text", x = threslabpos, y = min_err + min_err.sd/2,
                     label = "SD", hjust = -.25) +
            geom_vline(xintercept = res_phen$kopt, colour = "red", lty = "dashed") +
            annotate(geom = "label", x = res_phen$kopt, y = max(kopt_final$error),
                     label = paste("Optimal k =", res_phen$kopt), size = 6,
                     vjust = 1, hjust = -0.1, fontface = "bold") +
            labs(x = "Possible partitions", y = "Error rate") +
            theme(axis.title = element_text(size = 16),
                  axis.text = element_text(size = 12))
    }

ggsave("../plots/kopt_phen.png", kopt_phen)

##------------------------------
## Hierarchical clustering tree
##------------------------------

dend <- cov_phen %>%
    `class<-`("hclust") %>%
    ggtree(layout = "circular", branch.length = "none") %<+%
    unique(premat_phen %>% transmute(label = grpid, trait = str_trunc(trait, 55))) +
    geom_tiplab(aes(label = trait), size = 1.5) +
    xlim(0,25)

clades <- res_phen$ptree$cluster %>%
    data.frame %>%
    rownames_to_column %>%
    `names<-`(c("label", "Cluster")) %>%
    inner_join(dend$data) %>%
    mutate(Cluster = factor(Cluster)) %>%
    group_by(Cluster) %>%
    summarise(clade_num = MRCA(dend, label)) %>%
    mutate(clade_col = pals::cols25(nrow(.)))

gdend <- groupClade(dend, clades$clade_num) + aes(color = group) +
    scale_color_manual(values = c("#000000", clades$clade_col),
                       labels = c("", paste0("Cluster", 1:res_phen$kopt)),
                       guide = guide_legend(title = NULL, nrow = 1,
                                            keywidth = c(0, rep(1,res_phen$kopt)))) +
    theme(legend.position = "top")

ggsave("../plots/dend_phen.png", gdend)

##---------------------------------
## Cluster importance (color coded)
##---------------------------------

clusimp_phen <- with(res_phen$vsurf_ptree,
                     data.frame(Cluster = imp.mean.dec.ind,
                                clusnm = 1:length(imp.mean.dec.ind),
                                importance = imp.mean.dec,
                                importancesd = imp.sd.dec)) %>%
    inner_join(mutate(clades, Cluster = as.numeric(Cluster)))

clusimp_phen_plot <- clusimp_phen %>%
    mutate(Cluster = paste0("Cluster", Cluster)) %>%
    ggplot(aes(importance, reorder(Cluster, desc(clusnm)))) +
    geom_col(aes(fill = clade_col), color = "black", width = .5) +
    scale_fill_identity() +
    geom_errorbar(aes(xmin = importance - importancesd,
                      xmax = importance + importancesd),
                  width = 0.2, alpha = 0.5) +
    geom_hline(aes(yintercept = length(res_phen$vsel) + .5),
               linetype = "dashed") +
    geom_text(aes(x = max(importance), y = length(res_phen$vsel) + .5,
                  label = "Final model"), hjust = 1, vjust = -1, size = 3) +
    labs(title = "Importance score", x = "Average error rate decrease", y = NULL) +
    scale_x_continuous(limits = c(0, 
                                  with(clusimp_phen, 
                                       max(importance + importancesd) + 
                                           max(importance + importancesd) * .01)),
                       expand = c(0,0)) +
    theme(axis.title.x = element_text(size = 6),
          plot.title = element_text(size = 8, face = "bold"))

ggsave("../plots/clusimp_phen.png", clusimp_phen_plot)

##--------------------
## Cluster membership 
##--------------------

resvars_phen <- res_phen$ptree$cluster %>%
    data.frame %>%
    rownames_to_column() %>%
    `names<-`(c("grpid", "cluster")) %>%
    group_by(cluster) %>%
    group_modify(~{
        ## Extracting squared loading and correlation of traits
        dat <- res_phen$ptree$var[.y$cluster][[1]]
        if(is.null(dim(dat))){
            df <- data.frame(grpid = .x$grpid,
                             squared.loading = 1,
                             correlation = 1)
        } else {
            df <- dat %>%
                data.frame %>%
                rownames_to_column("grpid")
        }
        inner_join(.x, df, by = "grpid")
    }) %>%
    ungroup %>%
    inner_join(with(res_phen$vsurf_ptree,
                    data.frame(cluster = imp.mean.dec.ind,
                               importance = imp.mean.dec))) %>%
    mutate(inmod = cluster %in% res_phen$vsel) %>%
    inner_join(premat_phen) %>%
    arrange(match(cluster, res_phen$vsurf_ptree$imp.mean.dec.ind),
            desc(squared.loading))

rio::export(resvars_phen, "~/dva/files/resvars_phen.tsv")

## Summarized version
resvars_phen_sum <- resvars_phen %>%
    transmute(cluster, importance, included = inmod,
              squared.loading, ieu_id = id, trait,
              category, year, consortium, author, pmid, sample_size, nsnp,
              subtype = ifelse(disc == 1, "Discordant", "Concordant"),
              b_meta, se_meta, conf_low_meta, conf_high_meta,
              p_meta, Qv_meta = Qv, Qdf_meta = Qdf , Qp_meta = Qp, Tau_meta = Tau,
              D = b_comp, D_se = b_comp_se, D_p = p_comp) %>%
    unique()

rio::export(resvars_phen_sum, "~/dva/files/resvars_phen_sum.tsv")

## Quick check
resvars_phen_sum %>%
    filter(included, squared.loading > .5) %>%
    transmute(cluster, trait = str_trunc(trait, 50)) %>%
    data.frame

## Table for supplementary appendix
resvars_phen_sum %>%
    rio::export("~/dva/files/clusters.xlsx")

##---------------------------------------------
## Importance + structure of clusters selected
##---------------------------------------------

## Manual selection of main traits (up to SL > 0.25)
nms <- resvars_phen %>%
    filter(inmod, squared.loading > 0.25) %>%
    transmute(cluster, grpid, trait) %>%
    group_by(cluster) %>%
    group_map(~unique(paste(.x$grpid, .x$trait)))

nms

toselect <- paste0("v", c(2549, 2575, 1463, 1437,
                          1195, 136, 1355, 3128, 25,
                          3161, 3171, 3147, 1126, 1366, 3126, 970, 3094, 131, 3132,
                          3115, 3098, 3095,
                          219, 149, 1185,
                          159, 989, 2445, 223,
                          1212, 1213,
                          1977, 2341, 2145, 1938))

resvars_phen %>%
    select(grpid, trait) %>%
    unique %>%
    filter(grpid %in% toselect) %>%
    arrange(match(grpid, toselect)) %>%
    pull(trait)

## Editing names
newnms <- tibble(
    grpid = toselect,
    short_names = c("Arm fat mass", "Leg fat mass", "Body size age 10",
                    "Age at menarche", "HDL", "CHD",
                    "N meds taken", "ApoA", "Parental longevity", "SHBG", "Urate",
                    "GGT", "WHR", "Alcohol intake", "ALT", "AdipoQ", "#Leuco",
                    "ISI", "AST", "#Retic", "MCV", "#RBC", "Heel BMD",
                    "Total BMD", "Femoral neck BMD",
                    "Stroke", "Intracranial volume", "Extertional chest pain",
                    "FVC", "SBP", "DBP", "Total lean mass", "Arm lean mass",
                    "Basal met rate", "Leg lean mass")
)

## Data of traits selected
sigdiff <- resvars_phen %>%
    inner_join(newnms)

## Importance
imp_plot <- sigdiff %>%
    transmute(cluster, importance,
              short_names, squared.loading) %>%
    unique %>%
    inner_join(mutate(clades, cluster = as.numeric(Cluster))) %>%
    ggplot(aes(importance, short_names)) +
    geom_col(aes(fill = clade_col), width = 1) +
    scale_x_continuous(expand = c(0,0)) +
    scale_fill_identity() +
    ggforce::facet_col(vars(reorder(cluster, desc(importance))),
                       scales = "free_y", space = "free",
                       strip.position = "left") +
    labs(x = "Average error rate\ndecrease") +
    theme_classic() +
    theme(strip.text.y.left = element_text(angle = 0, color = "black", face = "bold"),
          strip.background = element_rect(color = "black", fill = "grey"),
          legend.position = "none",
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank())

## Squared loading
sl_plot <- sigdiff %>%
    select(cluster, importance, short_names, squared.loading) %>%
    unique %>%
    ggplot(aes(squared.loading, reorder(short_names, squared.loading))) +
    geom_col(width = .3) +
    scale_x_continuous(expand = c(0,0)) +
    ggforce::facet_col(vars(reorder(cluster, desc(importance))),
                       scales = "free_y", space = "free") +
    theme_linedraw() +
    theme(strip.text = element_blank(),
          strip.background = element_blank(),
          axis.title.y = element_blank(),
          plot.margin = unit(c(0,0,0,0), "cm")) +
    labs(x = "Squared loading")

## Pooled concordant and discordant effects
re_meta_plot <- sigdiff %>%
    transmute(cluster = paste("Cluster", cluster),
              importance, grpid, short_names, squared.loading,
              disc = factor(disc), b_meta, conf_low_meta, conf_high_meta) %>%
    ggplot(aes(b_meta, reorder(short_names, squared.loading), group = disc)) +
    geom_errorbar(aes(xmin = conf_low_meta, xmax = conf_high_meta),
                  alpha = .4, size = .4, width = .2,
                  position = position_dodge(width = .5)) +
    geom_point(aes(fill = disc), color = "black",
               alpha = .5, shape = 21,
               position = position_dodge(width = .5)) +
    scale_fill_discrete(labels = c("Concordant", "Discordant")) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    ggforce::facet_col(vars(reorder(cluster, desc(importance))),
                       scales = "free_y", space = "free") +
    theme_linedraw() +
    theme(strip.text = element_blank(),
          strip.background = element_blank(),
          legend.position = "top", legend.title = element_blank(),
          axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    labs(x = "Average per allele \u03b2", y = NULL)

## Final plot
clus_str_phen <- wrap_plots(sl_plot, re_meta_plot, imp_plot,
                            nrow = 1, widths = c(.25,.5,.25)) &
    theme(axis.title.x = element_text(size = 8),
          axis.text.x = element_text(size = 6))

ggsave("../plots/clus_str_phen.png")

##----------------
## Main variables
##----------------

## Main variables
mainvars_plot <- resvars_phen %>% 
    filter(grpid %in% c("v1212", "v1195", "v3161")) %>% 
    transmute(ref_rsid, disc = ifelse(disc == 1, "Discordant", "Concordant"), 
              trait = recode(trait, 
                             `systolic blood pressure` = "SBP",
                             `HDL cholesterol` = "HDL"), 
              zscore) %>% 
    pivot_wider(names_from = trait, values_from = zscore) %>%
    mutate(SHBG = cut(SHBG, c(-Inf, seq(-10,10,2.5), Inf), labels = FALSE), 
           alpha = scales::rescale(SHBG, c(0.4, 1))) %>%
    ggplot(aes(x = SBP, y = HDL, fill = SHBG)) +
    geom_point(aes(fill = disc, size = SHBG, 
                   alpha = alpha), color = "black", shape = 21) + 
    scale_alpha_identity() +
    scale_size_continuous(guide = guide_legend(title = "-\n\u2191\nSHBG\n\u2193\n+",
                                               label = FALSE,
                                               title.position = "right")) +
    geom_hline(yintercept = 0, linetype = "dashed") + 
        geom_vline(xintercept = 0, linetype = "dashed") +
        labs(x = "SBP", y = "HDL cholesterol", size = "SHBG", fill = NULL)

ggsave("../plots/mainvars_phen.png", mainvars_plot)

##--------------------
## Random Forest model
##--------------------

## Retraining RF multiple times
retrainrf <- lapply(1:1000,
                    function(x)
                        randomForest::randomForest(x = res_phen$ptree$scores,
                                                   y = res_phen$y,
                                                   ntree = 1000, proximity = TRUE))

## Extracting votes
all_votes <- retrainrf %>%
    lapply(function(x) data.frame(ref_rsid = mat_phen$ref_rsid,
                                  disc = mat_phen$disc,
                                  votes = x$votes[,2])) %>%
    bind_rows

## Extracting proximities
all_prox <- retrainrf %>%
    lapply(function(x) x$proximity) %>%
    {
        reduce(., `+`) / length(.)
    }

## ROC curve
rocres <- with(all_votes, roc(disc, votes))

rocplot <- ggroc(rocres) +
    coord_equal() +
    geom_abline(slope=1, intercept = 1, linetype = "dashed") +
    labs(x = "Specificity", y = "Sensitivity") +
    theme_minimal() +
    annotate("text", x = 0.25, y = 0.25,
             label = paste0("AUC=", round(rocres$auc, 2)))

ggsave("../plots/rocres.png")

## Selecting best threshold with ROC curve
best_thres <- coords(rocres, "best", transpose = TRUE)

## Random Forest classification
rfclass <- all_votes %>%
    group_by(ref_rsid) %>%
    summarise(rfprob = mean(votes),
              rfprob_sd = sd(votes)) %>%
    mutate(minvote = rfprob - rfprob_sd,
           maxvote = rfprob + rfprob_sd,
           rf_res = case_when("&"(minvote > best_thres[1],
                                  maxvote > best_thres[1]) ~ 1,
                              "&"(minvote < best_thres[1],
                                  maxvote < best_thres[1]) ~ 0,
                              TRUE ~ 0.5))

## SNP clustering
snpclus <- lapply(
    c(0,1),
    function(x){
        ## Proximity matrix
        snps <- mat_phen$ref_rsid[mat_phen$disc == x]
        dmat <- (1 -  all_prox[mat_phen$disc == x, mat_phen$disc == x])
        rownames(dmat) <- snps; colnames(dmat) <- snps
        ## Hierarchical clustering using proximity as distance
        d <- as.dist(dmat)
        h <- hclust(d, method = "ward.D2") ## Minimizing total within cluster-variance
        ## Optimal number of clusters using Dunn index
        nc_test <- NbClust(diss = d, distance = NULL, max.nc = length(snps)/2,
                           method = "ward.D2", index = "dunn")
        nc <- nc_test$Best.nc[1]
        nc_dat <- data.frame(nc_test$Best.partition) %>%
            rownames_to_column %>%
            `names<-`(c("label", "clus_nm"))
        ## Assessing cluster stability using bootstrapping (100 runs)
        stab_test <- clusterboot(d, clustermethod = disthclustCBI,
                                 k = nc, method = "ward.D2")
        stab_dat <- data.frame(clus_nm = 1:nc,
                               jaccard = round(stab_test$bootmean, 1))
        snpclusdat <- inner_join(nc_dat, stab_dat)
        return(list(snpclusdat = snpclusdat, hclusdat = h, nc = nc))
    }
) %>%
    `names<-`(c(0,1))

snpclus %>%
    map(~ .[["snpclusdat"]]) %>%
    bind_rows(.id = "disc") %>%
    rio::export("~/dva/files/snpclus.tsv")

## Heatmap with clustering based on proximity matrix
heatmap <- sigdiff %>%
    inner_join(rfclass) %>%
    mutate(disc = factor(disc, levels = c(1,0))) %>%
    group_by(disc) %>%
    group_map(~{
        ## SNP clusters
        snpclusinfo <- snpclus[[as.character(.y$disc)]]
        h <- snpclusinfo$hclusdat
        snp_c <- snpclusinfo$snpclusdat
        ## Joining cluster validation data
        dend_dat <- ggdendro::dendro_data(h) %>%
            "$"(labels) %>%
            inner_join(snp_c)
        ## Data for the heatmap of SNP effects
        hm_dat <- inner_join(.x, dend_dat, by = c("ref_rsid" = "label"))
        ## Colors for SNP clusters
        nc <- max(snp_c$clus_nm)
        snpcols <- unname(pals::alphabet2()[runif(nc, 1, 26)])
        ## Colors for trat clusters
        cluscol <- transmute(.x, Cluster = factor(cluster),
                             importance) %>%
            unique %>%
            inner_join(clades) %>%
            arrange(importance)
        ## Annotate Jaccard index
        jaccard_annot <- dend_dat %>%
            group_by(clus_nm) %>%
            filter(row_number() == floor(median(1:n())))
        ## Dendrogram
        hm_dend <- hm_dat %>%
            ggplot(aes(ref_rsid, "SNP clusters",
                       fill = factor(clus_nm))) +
            geom_raster(alpha = .75) +
            geom_text(data = jaccard_annot,
                      aes(x = label, y = "SNP clusters",
                          label = jaccard),
                      size = 2, hjust = .25) +
            scale_fill_manual(values = snpcols, guide = NULL) +
            scale_x_dendrogram(hclust = h, position = "top") +
            theme(axis.text.x = element_blank(),
                  axis.title = element_blank(),
                  axis.text.y = element_text(face = "bold.italic", size = 7))
        ## Heatmap
        hm_z <- hm_dat %>%
            mutate(ref_rsid = factor(ref_rsid, levels = dend_dat$label)) %>%
            ggplot(aes(ref_rsid,
                       interaction(reorder(short_names, squared.loading),
                                   reorder(cluster, importance)),
                       fill = zscore)) +
            geom_raster() +
            guides(y = "axis_nested") +
            scale_fill_gradient2(low = scales::muted("blue"),
                                 high = scales::muted("red"),
                                 guide = guide_colourbar(title = "- \u2190 Z-score \u2192 +",
                                                         title.position = "top",
                                                         title.theme = element_text(size = 8),
                                                         label = FALSE,
                                                         barheight = unit(5, "pt"))) +
            theme(ggh4x.axis.nestline.y = element_line(color = cluscol$clade_col,
                                                       size = 2.5),
                  title = element_blank(),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.text.y = element_text(size = 7),
                  axis.title = element_blank(),
                  legend.title.align = 0.5)
        if(.y$disc == 1){
            hm_z <- hm_z + guides(fill = NULL)
        } else {
            hm_dend <- hm_dend + theme(axis.text.y = element_blank(),
                                       axis.ticks.y = element_blank())
            hm_z <- hm_z + theme(axis.text.y = element_blank(),
                                 axis.ticks.y = element_blank(),
                                 ggh4x.axis.nestline.y = element_blank(),
                                 ggh4x.axis.nesttext.y = element_blank())
        }
        ## Random forest classification
        rf <- select(.x, ref_rsid, nearest_gene,
                        rfprob, rf_res) %>%
            unique %>%
            mutate(ref_rsid = factor(ref_rsid, levels = dend_dat$label),
                   snpgene = paste(nearest_gene, ref_rsid, sep = "-")) %>%
            ggplot(aes(reorder(snpgene, as.numeric(ref_rsid)),
                       "RF probability\n(discordant)", fill = rfprob)) +
            geom_raster() +
            geom_text(aes(label = ifelse(rf_res != .y$disc, "*", NA))) +
            scale_fill_gradient(high = "magenta",
                                low = "white", guide = NULL) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
                  axis.text.y = element_text(face = "bold.italic", size = 7),
                  axis.title = element_blank())
        if(.y$disc == 0){
            rf <- rf + theme(axis.text.y = element_blank(),
                             axis.ticks.y = element_blank())
        }
        ## Combining plots
        if(.y$disc == 0){
            p <- wrap_plots(hm_dend, hm_z, rf, ncol = 1, heights = c(.025, .95,.025),
                            guides = "collect") & theme(legend.position = "top")
        } else {
            p <- wrap_plots(hm_dend, hm_z, rf, ncol = 1, heights = c(.025, .95,.025)) &
                theme(legend.position = "none")
        }
        return(p)
    }) %>%
    wrap_plots(nrow = 1, widths = rev(table(res_phen$y))) &
    scale_y_discrete(expand = c(0,0)) &
    theme(plot.margin = unit(c(0,0,0,0), "cm"),
          panel.border = element_rect(colour = "black", fill = NA))

ggsave("../plots/phen_hm.png", heatmap)

## First two dimensions of proximity matrix
cmdscale(1 - all_prox, eig = TRUE) %>%
    "$"(points) %>%
    `colnames<-`(c("Dim1", "Dim2")) %>%
    data.frame %>%
    mutate(ref_rsid = mat_phen$ref_rsid) %>%
    inner_join(rfclass) %>%
    inner_join(mix, by = c("ref_rsid" = "rsid")) %>%
    mutate(labeltext = ifelse(rf_res != 0.5 & rf_res != disc,
                              paste(ref_rsid, nearest_gene, sep = "\n"), NA)) %>%
    ggplot(aes(Dim1, Dim2)) +
    geom_point(aes(color = factor(disc))) +
    scale_color_discrete(labels = c("Concordant", "Discordant")) +
    labs(color = NULL) +
    ggrepel::geom_text_repel(aes(label = labeltext), min.segment.length = 0,
                             box.padding = 1.5)

ggsave("../plots/mds_phen.png")
