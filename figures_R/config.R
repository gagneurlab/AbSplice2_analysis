# -*- coding: utf-8 -*-
suppressPackageStartupMessages({
  library(data.table)
  library(magrittr)
  library(ggplot2)
  library(cowplot)
  library(ggthemes)
  library(ggpubr)
  library(ggforce)
  library(ggrepel)
  library(gplots)
  library(png)
  library(ggvenn)
  library(eulerr)
  library(ggnewscale)
  require(ggbeeswarm)
  library(plyr)
  library(dplyr)
  library(RColorBrewer)
  library(heatmaply)
  library(ggplotify)
  library(UpSetR)
  library(Rmisc)
  library(ggupset)
  library(tidyverse, warn.conflicts = FALSE)
  library(stringr)
  library(tidyr)
})

fontsize <- 14
linesize <- 1
dpi=450

DATA_DIR <- "../workflow/data/results/source_data" # correct path
# DATA_DIR <- "../../AbSplice_analysis/workflow/data/results/source_data" # TODO: comment out

pr_curve_FRASER1_and_FRASER2_ground_truth <- file.path(DATA_DIR, "pr_curve_FRASER1_and_FRASER2_ground_truth.csv")
pr_curve_FRASER1_and_FRASER2_ground_truth_boxplot <- file.path(DATA_DIR, "pr_curve_FRASER1_and_FRASER2_ground_truth_boxplot.csv")

pr_curve_signed_features_SpliceAI_Pangolin <- file.path(DATA_DIR, "pr_curve_signed_features_SpliceAI_Pangolin.csv")
pr_curve_continous_median_n_impact <- file.path(DATA_DIR, "pr_curve_continous_median_n_impact.csv")
pr_curve_no_mmsplice_does_not_decrease_performance <- file.path(DATA_DIR, "pr_curve_no_mmsplice_does_not_decrease_performance.csv")
pr_curve_pangolin_affected_sites <- file.path(DATA_DIR, "pr_curve_pangolin_affected_sites.csv")

pr_curve_across_gene_expression_values <- file.path(DATA_DIR, "pr_curve_across_gene_expression_values.csv")
pr_curve_across_gene_expression_values_stats <- file.path(DATA_DIR, "pr_curve_across_gene_expression_values_stats.csv")

pr_curve_final_absplice2_model <- file.path(DATA_DIR, "pr_curve_final_absplice2_model.csv")
pr_curve_max_aggregation_decreases_performance <- file.path(DATA_DIR, "pr_curve_max_aggregation_decreases_performance.csv")

feature_importance_absplice1_fraser1_fraser2 <- file.path(DATA_DIR, "feature_importance_absplice1_fraser1_fraser2.csv")

absplice_model_colors = c(
  "AbSplice-DNA" = "#d62728",
  # pr_curve_FRASER1_and_FRASER2_ground_truth
  "AbSplice-DNA (trained on FRASER1, evaluated on FRASER1)" = 'red',
  'AbSplice-DNA (trained on FRASER2, evaluated on FRASER1)' = 'orange',
  'AbSplice-DNA (trained on FRASER1, evaluated on FRASER2)' = 'blue',
  'AbSplice-DNA (trained on FRASER2, evaluated on FRASER2)' = 'lightblue',
  
  # signed features
  "Pangolin Score" = 'grey',
  "SpliceAI Delta Score" = 'black',
  
  "AbSplice-DNA 1 (Pangolin max score)" = 'orange',
  "AbSplice-DNA 1 (Pangolin max score, signed)" = 'purple',
  "AbSplice-DNA 1 (Pangolin all scores, signed)" = 'red',
  
  "AbSplice-DNA 1 (SpliceAI max score)" = 'orange',
  "AbSplice-DNA 1 (SpliceAI max score, signed)"= 'purple',
  "AbSplice-DNA 1 (SpliceAI all scores, signed)" = 'red',
  
  # performance across gene expression values
  "AbSplice-DNA 1"= 'red',
  "AbSplice-DNA 2"= "#d62728",   
  "AbSplice-DNA with gene TPM" = 'orange', 
  "Pangolin"= 'blue',     
  "SpliceAI"= 'royalblue'
)



plot_pr_curve <- function(df, color_list, ylim=0.2, breaks_x=0.2, minor_breaks_x=0.2, breaks_y=0.1, minor_breaks_y=0.02, title='All tissues') {
  g <- (
    ggplot(df, aes(x=recall, y=precision, color=model))
    + geom_step(direction = "vh")
    + scale_color_manual(values=color_list, guide=guide_legend(reverse = TRUE))
    + theme_cowplot(font_size = fontsize)
    #+ background_grid(major = "xy", minor = "xy")
    + scale_x_continuous(limits = c(-0.01, 1.05), breaks = seq(0, 1, breaks_x), minor_breaks = seq(0, 1, minor_breaks_x))
    + scale_y_continuous(limits = c(-0.01, ylim), breaks = seq(0, 1, breaks_y), minor_breaks = seq(0, 1, minor_breaks_y))
    # + theme(
    #   legend.position = c(0.95, 0.95),
    #   legend.justification = c("right", "top")
    # )
    + labs(x='Recall', y='Precision', color='', title=title)
    + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    + theme(
      legend.position = "bottom",
      #legend.justification = "center"
    )
  )
  return(g)
}


plot_box_plot <- function(
    df, color_list, title="Across tissues",
    comparisons = list(
      c("AbSplice-DNA", "MMSplice + SpliceMap + Ψ_ref")
    ),
    coord_flip=FALSE,
    jitter=FALSE,
    x_lim=max(df$`Average Precision Score`) * 1.3
) {
  g <- (
    ggplot(df, aes(x=`Average Precision Score`, y=model, fill=model))
    + scale_fill_manual(values=color_list)
    + geom_boxplot()
    + xlim(0, x_lim)
    + scale_y_discrete(labels = function(x) {
      x %>% str_replace_all(" \\+", "\n+") %>% str_replace_all(" \\(", "\n(")
    })
    + theme_cowplot(font_size = fontsize)
    # + background_grid()
    + theme(legend.title = element_blank())
    + theme(legend.position = "right")
    + labs(
      x="auPRC",
      # x="area under the precision-recall curve",
      y=element_blank(),
      title=title
    )
    + {if (coord_flip == TRUE)coord_flip()}
    + stat_compare_means(
      comparisons = comparisons,
      method = "wilcox.test",
      paired = TRUE,
      method.args = list(alternative = "greater")
    )
    + {if (jitter == TRUE)geom_jitter(color="black", size=0.4, alpha=0.9)}
    #+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
    + theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    #axis.ticks.x = element_blank()
  )
  )
  return(g)
}

set_start_precision_zero <- function(df, chosen_models) {
  for (model_temp in chosen_models) {
    df_model = df[df[,df$model==model_temp]]
    df$precision[df$model==model_temp][nrow(df_model)] <- 0
  }
  return(df)
}

process_models_pr_curve <- function(df, chosen_models) {
  df <- df[df$model %in% chosen_models]
  #df$model <- factor(df$model, levels = chosen_models,ordered = TRUE)
  df <- set_start_precision_zero(df, chosen_models)
  return(df)
}

process_models_box_plot <- function(df, chosen_models) {
  df <- df[df$model %in% chosen_models]
  df$model <- factor(df$model, levels = chosen_models,ordered = TRUE)
  return(df)
}

tissue_mapping = c(
  'Adipose_Subcutaneous' = "Adipose Subcutaneous",
  'Adipose_Visceral_Omentum' = "Adipose Visceral Omentum",
  'Adrenal_Gland' = "Adrenal Gland",
  'Artery_Aorta' ="Artery Aorta",
  'Artery_Coronary' = "Artery Coronary",
  'Artery_Tibial' = "Artery Tibial",
  'Brain_Amygdala' = 'Brain Amygdala',
  'Brain_Anterior_cingulate_cortex_BA24' = 'Brain Ant. cing. cortex BA24',
  'Brain_Caudate_basal_ganglia' = 'Brain Caudate basal ganglia',
  'Brain_Cerebellar_Hemisphere' = 'Brain Cerebellar Hemisphere',
  'Brain_Cerebellum' = 'Brain Cerebellum',
  'Brain_Cortex' = 'Brain Cortex',
  'Brain_Frontal_Cortex_BA9' = 'Brain Frontal Cortex BA9',
  'Brain_Hippocampus' = 'Brain Hippocampus',
  'Brain_Hypothalamus' = 'Brain Hypothalamus',
  'Brain_Nucleus_accumbens_basal_ganglia' = 'Brain Nuc. basal ganglia',
  'Brain_Putamen_basal_ganglia' = 'Brain Putamen basal ganglia',
  'Brain_Spinal_cord_cervical_c_1' = 'Brain Spinal cord cervical c1',
  'Brain_Substantia_nigra' = 'Brain Substantia nigra',
  'Breast_Mammary_Tissue' = "Breast",
  'Colon_Sigmoid' = 'Colon Sigmoid',
  'Colon_Transverse' = 'Colon Transverse',
  'Esophagus_Gastroesophageal_Junction' = 'Esophagus GJ',
  'Esophagus_Mucosa' = 'Esophagus Mucosa',
  'Esophagus_Muscularis' = 'Esophagus Muscularis',
  'Heart_Atrial_Appendage' = 'Heart Atrial Appendage',
  'Heart_Left_Ventricle' = 'Heart Left Ventricle',
  'Kidney_Cortex' = 'Kidney Cortex',
  'Liver' = 'Liver',
  'Lung' = 'Lung',
  'Minor_Salivary_Gland' = 'Minor Salivary Gland',
  'Muscle_Skeletal' = 'Muscle Skeletal',
  'Nerve_Tibial' = 'Nerve Tibial',
  'Ovary' = 'Ovary',
  'Pancreas' = 'Pancreas',
  'Pituitary' = 'Pituitary',
  'Prostate' = 'Prostate',
  'Skin_Not_Sun_Exposed_Suprapubic' = 'Skin No Sun',
  'Skin_Sun_Exposed_Lower_leg' = 'Skin Sun',
  'Small_Intestine_Terminal_Ileum' = 'Small Intestine',
  'Spleen' = 'Spleen',
  'Stomach' = 'Stomach',
  'Testis' = 'Testis',
  'Thyroid' = 'Thyroid',
  'Uterus' = 'Uterus',
  'Vagina' = 'Vagina',
  'allBrainTissues' = 'Brain',
  'Cells_Cultured_fibroblasts' = 'Fibroblasts',
  'Cells_EBV_transformed_lymphocytes' = 'Lymphocytes',
  'Whole_Blood' = 'Blood'
)

