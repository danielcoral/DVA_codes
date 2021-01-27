## Plots of COVVSURF results

library(tidyverse)
library(patchwork)
library(cowplot)

res_phen <- rio::import("../files/res_phen.RData")

## Cluster importance
clusimp_phen <- data.frame(
    cluster = factor(res_phen$vsurf_ptree$varselect.thres),
    imp = res_phen$vsurf_ptree$imp.varselect.thres
) %>%
    filter(cluster %in% res_phen$vsurf_ptree$varselect.pred) %>%
    mutate(imp = 100 *(imp/sum(imp)),
           cluster = c("Hypertension", "Dyslipidaemia-CHD",
                       "Peripheral fat mass", "WHR"))

limits <- c(0, max(clusimp_phen$imp) + 5)

clusimp_phen <- clusimp_phen %>%
    ggplot() +
    geom_point(aes(imp, reorder(cluster, imp))) +
    geom_segment(aes(imp, reorder(cluster, imp),
                     xend = 0, yend = cluster)) +
    scale_x_continuous(limits = limits, expand = c(0,0)) +
    labs(x = "% CtA", y = "Clusters") +
    theme(panel.border = element_rect(colour = "black", fill = NA))

## ROC AUC
library(pROC)

rocres <- roc(res_phen$y, res_phen$rfsel$votes[,2])
roc_ci <- ci.se(rocres)
df_ci <- data.frame(x = as.numeric(rownames(roc_ci)),
                    lo = roc_ci[,1], up = roc_ci[,3])
rocplot <- ggroc(rocres) +
    geom_ribbon(data = df_ci, aes(x = x, ymin = lo, ymax = up),
                fill = "steelblue", alpha= 0.2) +
    coord_equal() +
    geom_abline(slope=1, intercept = 1, linetype = "dashed") +
    labs(x = "Specificity", y = "Sensitivity") +
    theme_minimal()

rocres

## Traits in the clusters
resvars_phen <- rio::import("../files/resvars_phen.tsv")

sigdiff <- resvars_phen %>%
    filter(squared.loading > 0.1, p_adj < 0.05) %>%
    ## Manual selection
    slice(3, 4, 11, 12, 13, 14, 17, 18, 27, 28, 33, 34,
          37, 38, 39, 40, 41, 42, 43, 44, 45, 46, ## First cluster
          47, 48, 51, 52, 53, 54, 57, 58, 59, 60, ## Second cluster
          61, 62, 63, 64, 67, 68, 69, 70,
          93, 94, 141, 142, 161, 162, 165, 166, ## Third cluster
          171, 172, 173, 174, 181, 182, ## Fourth cluster
          191, 192, 195, 196, 197, 198, 199, 200) %>%
    select(cluster, squared.loading, trait, disc, b, ci_lo, ci_up, p_adj)

## Shortening names
sigdiff$short_names <- rep(c(
    "Hypertension", "N meds taken", "N illnesses\n(non-cancer)", "Self-rated\nhealth",
    "SBP", "Alcohol intake\nfrequency", "DBP", "Adiponectin", "Exertional\nchest pain",
    "Insulin\ndisposition index", "Takes\nnitrate",
    "Cholesterol\nlowering meds", "Takes\nsimvastatin", "Takes\nlipitor",
    "Takes\naspirin", "HDL",
    "Weight", "Arm fat mass", "Leg fat mass", "BMI", "Basal metabolic\nrate",
    "Childhood\nobesity", "WHR\nmales", "Age at\nmenarche",
    "WHR adj\nsmoking", "WHR adj physical\nactivity", "WHR", "WHR\nfemales", 
    "Gout", "Daytime sleep", "Insulin sensitivity\nindex"),
    each = 2)

resphen_traits <- sigdiff %>%
    mutate(cluster = factor(recode(cluster,
                                   "1" = "Hypertension",
                                   "2" = "Dyslipidaemia",
                                   "4" = "Peripheral fat mass",
                                   "5" = "WHR"),
                            levels = c("Hypertension", "Dyslipidaemia",
                                       "Peripheral fat mass", "WHR")),
              disc = ifelse(disc == 1, "disc", "conc")) %>%
    group_by(cluster) %>%
    group_map(~list(
                  ## Importance of each trait
                  .x %>%
                  select(c(short_names, squared.loading)) %>%
                  ggplot(aes(squared.loading, "")) +
                  geom_point() +
                  geom_segment(aes(xend = 0, yend = "")) +
                  scale_x_continuous(limits = c(0,1), expand = c(0,0), n.breaks = 4) +
                  labs(x = "Squared Loading", title = .y$cluster) +
                  facet_wrap(reorder(short_names, desc(squared.loading)) ~ .,
                             ncol = 1, strip.position = "left") +
                  theme(axis.title.y = element_blank(),
                        axis.text.y = element_blank(),
                        axis.ticks.y = element_blank(),
                        axis.title.x = element_text(size = 6),
                        axis.text.x = element_text(size = 6),
                        strip.text.y.left = element_text(size = 7, angle = 0),
                        strip.background = element_blank(),
                        plot.title.position = "plot",
                        plot.title = element_text(size = 8, face = "bold",
                                                  hjust = 0.5)),
                  ## Effects of concordant/discordant SNPs for each trait
                  ggplot(.x, aes(b, disc)) +
                  geom_errorbar(aes(xmin = ci_lo, xmax = ci_up),
                                alpha = .4, size = .4, width = .2) +
                  geom_point(aes(colour = disc), alpha = .5) +
                  geom_point(shape = 1, col = "black", alpha = .4) +
                  scale_colour_manual(values = c("red", "blue"),
                                      labels = c("Concordant", "Discordant")) +
                  guides(colour = guide_legend(title = NULL)) +
                  geom_vline(xintercept = 0, linetype = "dashed") +
                  xlab("Pooled \u03b2") +
                  facet_wrap(reorder(short_names, desc(squared.loading)) ~ ., ncol = 1) +
                  theme(axis.text.y = element_blank(),
                        axis.ticks.y = element_blank(),
                        axis.title.y = element_blank(),
                        axis.text.x = element_text(size = 6),
                        axis.title.x = element_text(size = 6),
                        strip.text = element_blank(),
                        strip.background = element_blank())        
              ) %>%
                  wrap_plots(nrow = 1, widths = c(3,7))) %>%
    wrap_plots(ncol = 2) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

## Join plots
covsurf_phen <- plot_grid(
    plot_grid(clusimp_phen, rocplot, labels = c("A", "B")),
    resphen_traits, labels = c("", "C"), ncol = 1,
    rel_heights = c(1,3.5)
)

ggsave("../plots/covsurf_phen.png", covsurf_phen, height = 10)
