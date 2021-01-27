## Plots of CoVVSURF steps
library(tidyverse)

cov_phen <- rio::import("../files/cov_phen.RData")
res_phen <- rio::import("../files/res_phen.RData")                                               
id_phen <- rio::import("../files/id_phen.txt")

## Dendrogram
library(ggtree)

dend <- cov_phen %>%
    `class<-`("hclust") %>%
    ggtree(layout = "circular", branch.length = "none") %<+%
    select(id_phen, label = id, trait) +
    geom_tiplab(aes(label = trait), size = .8) +
    xlim(0,25)

clades <- res_phen$ptree$cluster %>%
    data.frame %>%
    rownames_to_column %>%
    `names<-`(c("label", "Cluster")) %>%
    inner_join(dend$data) %>%
    mutate(Cluster = factor(Cluster)) %>%
    group_by(Cluster) %>%
    summarise(clade_num = MRCA(dend, label)) %>%
    mutate(clade_col = RColorBrewer::brewer.pal(length(clade_num), "Set2"))

gdend <- groupClade(dend, clades$clade_num) + aes(color = group) +
    scale_color_manual(values = c("#000000", clades$clade_col)) +
    theme(legend.position = "none")

ggsave("../plots/dend_phen.png", gdend)

### Optimal number of clusters through RF error rate
err_df <- cbind(res_phen$kval, res_phen$oob) %>%
    data.frame %>%
    `names<-`(c("k", "error", "error_sd"))
min_err <- min(err_df$error)
min_err.sd <- err_df$error_sd[which.min(err_df$error)]
kopt_thresh <- min_err + min_err.sd
kopt_phen <- err_df %>%
    ggplot(aes(k, error)) + geom_line() + 
    geom_hline(yintercept = min_err, colour = "gray", lty = "dashed") +
    annotate("text", x = max(err_df$k), y = min_err, label = "Lowest error",
             vjust = 1.5, hjust = 1) +
    geom_hline(yintercept = kopt_thresh, colour = "maroon", lty = "dashed") +
    annotate("text", x = max(err_df$k), y = kopt_thresh, label = "Threshold",
             vjust = -0.5, hjust = 1) +
    geom_segment(aes(x = max(k), y = min_err, xend = max(k), yend = kopt_thresh),
                 arrow = arrow(ends = "both", length = unit(.2, "cm"))) +
    annotate("text", x = max(err_df$k), y = min_err + min_err.sd/2,
             label = "Lowest error SD", hjust = 1.1) +
    geom_vline(xintercept = res_phen$kopt, colour = "red", lty = "dashed") +
    annotate(geom = "label", x = res_phen$kopt, y = max(err_df$error),
             label = paste("Optimal k =", res_phen$kopt), size = 6,
             vjust = 1, hjust = -0.1, fontface = "bold") +
    labs(x = "Number of clusters", y = "Error rate") +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 12))

ggsave("../plots/kopt.png", kopt_phen)

## VSURF plots
library(patchwork)

vsurf_thres1 <- clades %>%
    mutate(measure = "VI mean",
           Value = res_phen$vsurf_ptree$imp.mean.dec) %>%
    ggplot(aes(Cluster, Value, group = 1)) +
    geom_line() + geom_point(aes(colour = clade_col)) +
    scale_colour_identity(guide = "none") + 
    geom_hline(yintercept = res_phen$vsurf_ptree$min.thres,
               linetype = "dashed", color = "red") +
    annotate("text", x = 1, y = res_phen$vsurf_ptree$min.thres,
             label = "Threshold 1", vjust = -0.5, hjust = -0.25, size = 3) +
    labs(y = "Error rate increase") +
    facet_wrap(~ measure)

vsurf_thres2 <- rbind(mutate(clades, measure = "VI SD",
                             Value = res_phen$vsurf_ptree$imp.sd.dec),
                      mutate(clades, measure = "Decision tree",
                             Value = res_phen$vsurf_ptree$pred.pruned.tree)) %>%
    mutate(measure = factor(measure, levels = unique(measure))) %>%
    ggplot(aes(Cluster, Value, group = 1)) +
    geom_line() + geom_point(aes(colour = clade_col)) +
    scale_colour_identity(guide = "none") +
    geom_hline(yintercept = res_phen$vsurf_ptree$min.thres,
               linetype = "dashed", color = "red") +
    labs(y = "Error rate increase - SD") +
    facet_wrap(~ measure)

vsurf_interp <- clades %>%
    filter(Cluster %in% res_phen$vsurf_ptree$varselect.thres) %>%
    mutate(measure = "Nested models",
           OOB_error = res_phen$vsurf_ptree$err.interp) %>%
    ggplot(aes(Cluster, OOB_error, group = 1)) +
    geom_line() + geom_point(aes(colour = clade_col)) +
    scale_colour_identity(guide = "none") +
    geom_vline(xintercept = length(res_phen$vsurf_ptree$varselect.interp),
               linetype = "dashed", color = "red") +
    geom_rect(aes(xmin = length(res_phen$vsurf_ptree$varselect.interp),
                  xmax = Inf, ymin = -Inf, ymax = Inf, fill = "darksalmon"),
              alpha = 0.05) +
    scale_fill_identity() +
    geom_hline(aes(yintercept = min(OOB_error)), colour = "gray", lty = "dashed") +
    annotate("text", x = 1, y = min(res_phen$vsurf_ptree$err.interp),
             label = "Lowest error", vjust = 1, hjust = 0.25, size = 2) +
    geom_hline(aes(yintercept = min(OOB_error) + res_phen$vsurf_ptree$sd.min),
               colour = "maroon", lty = "dashed") +
    annotate("text", x = 1,
             y = min(res_phen$vsurf_ptree$err.interp) + res_phen$vsurf_ptree$sd.min,
             label = "Threshold 2", vjust = -0.5, hjust = 0.25, size = 3) +
    annotate(geom = "label", x = length(res_phen$vsurf_ptree$varselect.interp),
             y = max(res_phen$vsurf_ptree$err.interp),
             label = "Optimal model", size = 3,
             vjust = 1, hjust = 1.1, fontface = "bold") +
    labs(y = "OOB error rate") +
    facet_wrap(~ measure)

vsurf_pred <- clades %>%
    filter(Cluster %in% res_phen$vsel) %>%
    mutate(measure = "Final model", OOB_error = res_phen$vsurf_ptree$err.pred) %>%
    mutate(diff = round(OOB_error - lag(OOB_error), 3)) %>%
    ggplot(aes(Cluster, OOB_error, group = 1)) +
    geom_line() + geom_point(aes(colour = clade_col)) +
    scale_colour_identity(guide = "none") +
    geom_text(aes(label = diff), hjust = 1.1, vjust = 1.1, size = 3) +
    annotate("text", x = 1, y = 0, vjust = -0.5, hjust = 0.25, size = 3,
             label = paste0("Threshold 3 = ",
                            round(res_phen$vsurf_ptree$mean.jump, 3))) +
    facet_wrap(~ measure) +
    theme(axis.title.y = element_blank())

vsurf <- wrap_plots(vsurf_thres1, vsurf_thres2,
                    vsurf_interp, vsurf_pred,
                    ncol = 2) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = 'bold'))

ggsave("../plots/vsurf.png", vsurf)
