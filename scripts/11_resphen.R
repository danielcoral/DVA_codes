## COVVSURF results - Phenome matrix

library(tidyverse)
library(ggtree)
library(patchwork)
library(ggh4x)
library(NbClust)
library(fpc)
library(randomForest)
library(grid)

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
                       guide = guide_legend(title = "Clusters"))

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
    labs(x = "Importance score", y = NULL) +
    scale_x_continuous(limits = c(0, 
                                  with(clusimp_phen, 
                                       max(importance + importancesd) + 
                                           max(importance + importancesd) * .01)),
                       expand = c(0,0)) +
    theme_light()

ggsave("../plots/clusimp_phen.png", clusimp_phen_plot)

## Nested random forests
with(res_phen$vsurf_ptree,
     data.frame(clusnm = 1:res_phen$kopt,
                clus = imp.mean.dec.ind,
                error_rate = c(err.interp,
                               rep(NaN, res_phen$kopt - length(err.interp)))) %>%
     mutate(mod = paste0(c("", rep("+",res_phen$kopt - 1)), "cluster", clus),
            minerr = min(error_rate, na.rm = TRUE), minerrsd = sd.min,
            selected = clus %in% res_phen$vsel,
            lastvar = last(mod[selected]))) %>%
    ggplot(aes(error_rate, reorder(mod, desc(clusnm)), group = 1)) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = unique(lastvar), ymax = Inf),
              fill = "lightgreen") +
    geom_hline(aes(yintercept = unique(lastvar)), linetype = "dashed", color = "darkgreen") +
    geom_line(orientation = "y") + geom_point() +
    geom_text(aes(max(error_rate, na.rm = TRUE), unique(lastvar)),
              label = "Final model", vjust = -1, hjust = 1) +
    labs(x = "Average error rate", y = "Nested models") +
    theme_light()

ggsave("../plots/nestedmodels.png")

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
                    "Stroke", "Intracranial volume", "Exertional chest pain",
                    "FVC", "SBP", "DBP", "Total lean mass", "Arm lean mass",
                    "Basal met rate", "Leg lean mass")
)

## Data of traits selected
sigdiff <- resvars_phen %>%
    inner_join(newnms)

##--------------------
## Random Forest model
##--------------------

## Obtaining RF probabilities assigned to each SNP
rf_res <- select(mat_phen, ref_rsid) %>%
    mutate(prob = predict(res_phen$rfsel, type = "prob")[,2],
           pred = predict(res_phen$rfsel, type = "response"))

## Obtaining proximities
proxmat <- predict(res_phen$rfsel, res_phen$ptree$scores[,res_phen$vsel], proximity = TRUE)$proximity

## SNP clustering
snpclus <- lapply(
    c(0,1),
    function(x){
        ## Proximity matrix
        snps <- mat_phen$ref_rsid[mat_phen$disc == x]
        nearest_gene <- mix$nearest_gene[match(snps, mix$rsid)]
        snps <- paste0(nearest_gene, " (", snps, ")")
        dmat <- (1 - proxmat[mat_phen$disc == x, mat_phen$disc == x])
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

## Importance of clusters selected
selclusimp <- sigdiff %>%
    transmute(cluster = factor(cluster, levels = unique(cluster)),
              grpid, importance) %>%
    unique()

colsizes <- selclusimp %>%
    group_by(cluster) %>%
    summarise(s = n())

selclusimp_p <- selclusimp %>%
    ggplot(aes(grpid, importance)) +
    geom_col(width = 1) +
    scale_y_reverse(expand = c(0,0), position = "right") +
    facet_wrap(~cluster, nrow = 1, scales = "free_x") +
    theme_light() +
    labs(y = "Importance\nscore") +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 4),
          axis.title.y.right = element_text(angle = 0, vjust = 0.5, size = 8),
          strip.background = element_blank(),
          strip.text = element_blank()) +
    force_panelsizes(cols = colsizes$s)

## Squared loading of traits
sq <- sigdiff %>%
    transmute(cluster = factor(cluster, levels = unique(cluster)),
              short_names = factor(short_names, levels = unique(short_names)),
              squared.loading) %>%
    unique %>%
    ggplot(aes(short_names, squared.loading)) +
    geom_col() +
    scale_y_reverse(expand = c(0,0), position = "right") +
    facet_wrap(~cluster, nrow = 1, scales = "free_x") +
    theme_light() +
    labs(y = "Squared\nloading") +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 4),
          axis.title.y.right = element_text(angle = 0, vjust = 0.5, size = 8),
          strip.background = element_blank(),
          strip.text = element_blank()) +
    force_panelsizes(cols = colsizes$s)

## Heatmap
rowsize <- sigdiff %>%
    transmute(disc = factor(disc, levels = c(1,0)), ref_rsid) %>%
    unique %>%
    group_by(disc) %>%
    summarise(n = n())

heatmap <- sigdiff %>%
    mutate(cluster = factor(cluster, levels = unique(cluster)),
           short_names = factor(short_names, levels = unique(short_names)),
           disc = factor(disc, levels = c(1,0))) %>%
    inner_join(rf_res) %>%
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
        hm_dat <- .x %>%
            mutate(ref_rsid = paste0(nearest_gene, " (", ref_rsid, ")")) %>%
            inner_join(dend_dat, by = c("ref_rsid" = "label"))
        ## Colors for SNP clusters
        nc <- max(snp_c$clus_nm)
        snpcols <- unname(pals::alphabet2()[runif(nc, 1, 26)])
        ## Annotate Jaccard index
        ji <- dend_dat %>%
            group_by(clus_nm) %>%
            filter(row_number() == floor(median(1:n())))
        return(list(disc = .y$disc, dat = .x, snpclusinfo = snpclusinfo, h = h, snp_c = snp_c,
                    dend_dat = dend_dat, hm_dat = hm_dat, nc = nc,
                    snpcols = snpcols, ji = ji))
    })

## SNP cluster dendrograms
snpdendp <- heatmap %>%
    map(~{
        dendp <- ggplot(.$hm_dat,
                        aes("JI", ref_rsid, fill = factor(clus_nm))) +
            geom_raster(alpha = .5) +
            geom_text(data = .$ji,
                      aes(x = "JI", y = label, label = jaccard), size = 2) +
            scale_fill_manual(values = .$snpcols, guide = NULL) +
            scale_y_dendrogram(hclust = .$h) +
            theme(axis.text.y = element_text(size = 5), axis.title.x  = element_blank(),
                  plot.margin = unit(c(0,0,0,0), "cm"),
                  panel.background = element_blank())
        if(.$disc == 1){
            dendp <- dendp + scale_x_discrete(expand = c(0,0), position = "top") +
                labs(y = "Discordant")
        } else {
            dendp <- dendp + scale_x_discrete(expand = c(0,0)) +
                labs(y = "Concordant") +
                theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
        }
        dendp
    }) %>%
    wrap_plots(ncol = 1, heights = rowsize$n)

## Random forest classification
rfp <- heatmap %>%
    map(~{
        hm_dat <- .$hm_dat
        dend_dat <- .$dend_dat
        disc <- .$disc
        p <- select(hm_dat, ref_rsid, prob, pred) %>%
            unique %>%
            mutate(ref_rsid = factor(ref_rsid, levels = dend_dat$label)) %>%
            ggplot(aes("RF", ref_rsid, fill = prob)) +
            geom_raster() +
            geom_text(aes(label = ifelse(pred != disc, "X", NA)), size = 2) +
            scale_fill_gradient(high = "magenta", low = "white", guide = NULL) +
            theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                  axis.title = element_blank(), plot.margin = unit(c(0,0,0,0), "cm"))
        if(disc == 1){
            p <- p + scale_x_discrete(expand = c(0,0), position = "top")
        } else {
            p <- p + scale_x_discrete(expand = c(0,0)) +
                theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
        }
        p
    }) %>%
    wrap_plots(ncol = 1, heights = rowsize$n)

## Heatmap
hmp <- heatmap %>%
    map(~{
        hm_dat <- .$hm_dat
        dend_dat <- .$dend_dat
        disc <- .$disc
        hm_z <- hm_dat %>%
            mutate(ref_rsid = factor(ref_rsid, levels = dend_dat$label)) %>%
            ggplot(aes(short_names, ref_rsid, fill = zscore)) +
            geom_raster() +
            scale_fill_gradient2(low = scales::muted("blue"),
                                 high = scales::muted("red"),
                                 guide = guide_colourbar(title = "- \u2190 Z-score \u2192 +",
                                                         title.position = "top",
                                                         title.theme = element_text(size = 8),
                                                         label = FALSE,
                                                         barheight = unit(5, "pt"))) +
            facet_wrap(~cluster, nrow = 1, scales = "free_x") +
            theme(axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.title = element_blank(),
                  strip.background = element_blank(),
                  strip.text = element_blank()) +
            force_panelsizes(cols = colsizes$s)
        if(disc == 0){
            hm_z <- hm_z + theme(axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(),
                                 legend.position = "none")
        } else {
            hm_z <- hm_z +
                scale_x_discrete(position = "top") +
                theme(axis.text.x = element_text(size = 6, angle = 45, hjust = .01),
                      axis.text.x.top = element_text(vjust = .01),
                      legend.position = "top", legend.title.align = 0.5)
        }              
        hm_z
    }) %>%
    wrap_plots(ncol = 1, heights = rowsize$n) &
    theme(plot.margin = unit(c(0,0,0,0), "cm"))

## Combining plots
heatmap_c <- wrap_plots(snpdendp, rfp, hmp,
                        plot_spacer(), plot_spacer(), sq,
                        plot_spacer(), plot_spacer(), selclusimp_p,
                        ncol = 3, widths = c(.04, .03, .93),
                        heights = c(.9, .1, .1)) &
    theme(panel.spacing = unit(0, "cm"),
          plot.margin = unit(c(0,0.025,0,0.025), "cm"),
          panel.border = element_rect(colour = "black", fill = NA))

ggsave("../plots/phen_hm.png", heatmap_c, width = 170, height = 200, units = "mm")


