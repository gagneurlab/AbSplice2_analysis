source("config.R")

# ===================== AbSplice2 overall ===================== 

df <- fread(pr_curve_final_absplice2_model)
df_rna <- fread(file.path(DATA_DIR, 'pr_curve_AbSplice_RNA.csv'))

df_rna$model[df_rna$model == "AbSplice2_RNA"] <- "AbSplice2-RNA"
df_rna <- df_rna[df_rna$model == "AbSplice2-RNA", ]

df <- rbind(df, df_rna)

colors_here <- c(
  "AbSplice-DNA 1 (trained on FRASER1, evaluated on FRASER1)" = "#1f77b4",  # Blue
  "AbSplice-DNA 1 (trained on FRASER1, evaluated on FRASER2)" = "#ff7f0e",  # Orange
  "AbSplice-DNA (SpliceAI -> Pangolin)" = "#2ca02c",  # Green
  "AbSplice-DNA 2" = "#d62728",  # Strong Red
  "AbSplice2-RNA" = "#9467bd",  # Purple
  "SpliceAI" = "#b3d9ff",
  "Pangolin" = "#ffd699" 
)
chosen_models <- c(
  "SpliceAI",
  "Pangolin",
  "AbSplice-DNA 1 (trained on FRASER1, evaluated on FRASER1)",
  "AbSplice-DNA 1 (trained on FRASER1, evaluated on FRASER2)",
  "AbSplice-DNA (SpliceAI -> Pangolin)",
  "AbSplice-DNA 2",
  "AbSplice2-RNA"
)

# Dummy label dictionary (can be edited later)
label_dict <- c(
  "SpliceAI" = "SpliceAI",
  "Pangolin" = "Pangolin",
  "AbSplice-DNA 1 (trained on FRASER1, evaluated on FRASER1)" = "AbSplice1 \n(trained on FRASER1, \nevaluated on FRASER1)",
  "AbSplice-DNA 1 (trained on FRASER1, evaluated on FRASER2)" = "AbSplice1 \n(trained on FRASER1, \nevaluated on FRASER2)",
  "AbSplice-DNA (SpliceAI -> Pangolin)" = "AbSplice (Pangolin)",
  "AbSplice-DNA 2" = "AbSplice2",
  "AbSplice2-RNA" = "AbSplice2-RNA"
)

df_here <- df[model %in% chosen_models]
df_here <- process_models_pr_curve(df_here, chosen_models)
df_here$model <- factor(df_here$model, levels = chosen_models)

# Remap color keys to match model_label values
color_list <- setNames(colors_here[names(label_dict)], label_dict[names(label_dict)])


df_here$model_label <- factor(label_dict[as.character(df_here$model)],
                              levels = label_dict[chosen_models])



plot_pr_curve_GTEx <- function(df, color_list, ylim=0.2, breaks_x=0.2, minor_breaks_x=0.2, breaks_y=0.1, minor_breaks_y=0.02, title='All tissues') {
  g <- (
    ggplot(df, aes(x=recall, y=precision, color=model_label))
    + geom_step(direction = "vh", linewidth=0.9)
    + scale_color_manual(values=color_list, guide=guide_legend(reverse = TRUE, ncol=1))
    + theme_cowplot(font_size = fontsize)
    + scale_x_continuous(limits = c(-0.01, 1.05), breaks = seq(0, 1, breaks_x), minor_breaks = seq(0, 1, minor_breaks_x))
    + scale_y_continuous(limits = c(-0.01, ylim), breaks = seq(0, 1, breaks_y), minor_breaks = seq(0, 1, minor_breaks_y))
    + theme(
      legend.position = c(0.95, 0.95),
      legend.justification = c("right", "top")
    )
    + labs(x='Recall', y='Precision', color='', title=title)
    + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  )
  return(g)
}

fig_pr_curve_gtex <- plot_pr_curve_GTEx(
  df_here,
  ylim=1.0,
  breaks_y=0.2,
  color_list,
  title='Aberrant splicing prediction \n(GTEx, all tissues)',
)
fig_pr_curve_gtex


# ========= boxplot

df_rna <- fread(file.path(DATA_DIR, 'pr_curve_AbSplice_RNA.csv'))
df <- fread(file.path(DATA_DIR, 'boxplot_all_tissues_GTEx.csv'))
df_rna <- fread(file.path(DATA_DIR, 'boxplot_all_tissues_GTEx_AbSplice_RNA.csv'))

df_rna$model[df_rna$model == "AbSplice2_RNA"] <- "AbSplice2-RNA"
df_rna <- df_rna[df_rna$model == "AbSplice2-RNA", ]
df <- rbind(df, df_rna)

df <- df[model %in% chosen_models]
df <- process_models_box_plot(df, unique(df$model))
df$model <- factor(df$model, levels = chosen_models)

color_list <- setNames(colors_here[names(label_dict)], label_dict[names(label_dict)])
df$model_label <- factor(label_dict[as.character(df$model)],
                              levels = label_dict[chosen_models])

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
    ggplot(df, aes(x=`Average Precision Score`, y=model, fill=model_label))
    + scale_fill_manual(values=color_list, guide=guide_legend(reverse = TRUE, nrow=4))
    + geom_boxplot()
    + xlim(0, x_lim)
    + theme_cowplot(font_size = fontsize)
    + theme(legend.title = element_blank())
    + theme(legend.position = "none")
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
      method.args = list(alternative = "greater", exact = FALSE)  # <-- Fix: Use approximation for ties
    )
    + {if (jitter == TRUE)geom_jitter(color="black", size=0.4, alpha=0.9)}
    + theme(
      axis.title.x = element_blank(),
      axis.text.x  = element_blank(),
    )
  )
  return(g)
}

fig_box <- plot_box_plot(
  df,
  color_list, 
  comparisons=list(
    # c("AbSplice-DNA 1 (trained on FRASER1, evaluated on FRASER1)", "AbSplice-DNA 1 (trained on FRASER1, evaluated on FRASER2)"),
    c("Pangolin", "SpliceAI"),
    # c("AbSplice-DNA (SpliceAI -> Pangolin)", "AbSplice-DNA 2"),
    c("AbSplice-DNA 2", "AbSplice-DNA 1 (trained on FRASER1, evaluated on FRASER1)"),
    c("AbSplice2-RNA", "AbSplice-DNA 2")
  ),
  coord_flip=TRUE,
  title='Aberrant splicing prediction \n(GTEx, across tissues)'
)
fig_box

## ============ GVEX
plot_pr_curve_facet <- function(df, color_list, ylim=0.2, breaks_x=0.2, minor_breaks_x=0.2, breaks_y=0.1, minor_breaks_y=0.02, title='All tissues') {
  g <- (
    ggplot(df, aes(x=recall, y=precision, color=model))
    + geom_step(direction = "vh", linewidth=0.9)
    + scale_color_manual(values=color_list, guide=guide_legend(reverse = TRUE, nrow=3))
    + theme_cowplot(font_size = fontsize)
    + facet_wrap('dataset_label')
    + scale_x_continuous(limits = c(-0.01, 1.05), breaks = seq(0, 1, breaks_x), minor_breaks = seq(0, 1, minor_breaks_x))
    + scale_y_continuous(limits = c(-0.01, ylim), breaks = seq(0, 1, breaks_y), minor_breaks = seq(0, 1, minor_breaks_y))
    + labs(x='Recall', y='Precision', color='', title=title)
    + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    + theme(legend.position = 'none')
  )
  return(g)
}

df <- fread(file.path(DATA_DIR, 'pr_curve_gvex_fraser2_Brain_Frontal_Cortex_BA9_compare_to_GTEx.csv'))

dataset_map <- c(
  "GTEx Brain_Frontal_Cortex_BA9" = "GTEx Brain Frontal Cortex",
  "GVEx Brain_Frontal_Cortex_BA9" = "GVEx Brain Frontal Cortex"
)

df <- df %>%
  mutate(dataset_label = dataset_map[dataset])

colors_here <- c(
  "AbSplice-DNA 1" = "#ff7f0e",  # Orange
  "AbSplice-DNA 2" = "#d62728",  # Strong Red
  "SpliceAI" = "#b3d9ff",
  "Pangolin" = "#ffd699" 
)
chosen_models <- c(
  "SpliceAI",
  "Pangolin",
  "AbSplice-DNA 1",
  "AbSplice-DNA 2"
)

df_here <- df[model %in% chosen_models]
df_here$model <- factor(df_here$model, levels = chosen_models)

color_list = colors_here[unique(df_here$model)]

fig_pr_curve_gvex <- plot_pr_curve_facet(
  df_here,
  ylim=1.0,
  breaks_y=0.2,
  color_list,
  title='Aberrant splicing prediction \n(different cohorts)',
)
fig_pr_curve_gvex


# ============== aberrant underexpression

plot_pr_curve <- function(df, color_list, ylim=0.2, breaks_x=0.2, minor_breaks_x=0.2, breaks_y=0.1, minor_breaks_y=0.02, title='All tissues') {
  g <- (
    ggplot(df, aes(x=recall, y=precision, color=model))
    + geom_step(direction = "vh", linewidth=0.5)
    + scale_color_manual(values=color_list, guide=guide_legend(reverse = TRUE, nrow=3))
    + theme_cowplot(font_size = fontsize)
    + scale_x_continuous(limits = c(-0.01, 1.05), breaks = seq(0, 1, breaks_x), minor_breaks = seq(0, 1, minor_breaks_x))
    + scale_y_continuous(limits = c(-0.01, ylim), breaks = seq(0, 1, breaks_y), minor_breaks = seq(0, 1, minor_breaks_y))
    + labs(x='Recall', y='Precision', color='', title=title)
    + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    + theme(legend.position = 'none')
  )
  return(g)
}

df <- fread(file.path(DATA_DIR, 'pr_curve_gtex_underexpression_outliers.csv'))

colors_here <- c(
  "AbSplice-DNA 1" = "#ff7f0e",  # Orange
  "AbSplice-DNA (SpliceAI -> Pangolin)" = "#2ca02c",  # Green
  "AbSplice-DNA 2" = "#d62728",  # Strong Red
  "SpliceAI" = "#b3d9ff",
  "Pangolin" = "#ffd699" 
)
chosen_models <- c(
  "SpliceAI",
  "Pangolin",
  "AbSplice-DNA 1",
  "AbSplice-DNA (SpliceAI -> Pangolin)",
  "AbSplice-DNA 2"
)

df_here <- df[model %in% chosen_models]
# df_here <- process_models_pr_curve(df_here, chosen_models)
df_here$model <- factor(df_here$model, levels = chosen_models)

color_list = colors_here[unique(df_here$model)]

# Load data
df <- fread(file.path(DATA_DIR, 'AbExp_performance.csv'))

# Define models and colors
absplice_models <- c("AbSplice1", "AbSplice2")
abexp_models <- c("AbExp \n(AbSplice1)", "AbExp \n(AbSplice2)")

colors_here <- c(
  "AbSplice1" = "#ff7f0e",
  "AbSplice2" = "#d62728",
  "AbExp \n(AbSplice1)" = "#dfb23c",
  "AbExp \n(AbSplice2)" = "#65274d"
)

# Filter separately
df_absplice <- df %>%
  filter(model_name %in% absplice_models) %>%
  mutate(model_name = factor(model_name, levels = absplice_models))

df_abexp <- df %>%
  filter(model_name %in% abexp_models) %>%
  mutate(model_name = factor(model_name, levels = abexp_models))

# Plot 1: AbSplice
plot_absplice <- ggplot(df_absplice, aes(y = model_name, x = auc, fill = model_name)) +
  geom_boxplot() +
  scale_fill_manual(values = colors_here) +
  stat_compare_means(
    comparisons = list(c("AbSplice2", "AbSplice1")),
    method = "wilcox.test",
    paired = TRUE,
    method.args = list(alternative = "greater", exact = FALSE)
  ) +
  labs(y = NULL, x = "auPRC", title = "") +
  theme_cowplot(font_size = fontsize) +
  theme(legend.position = "none")
plot_absplice


# Plot 2: AbExp
plot_abexp <- ggplot(df_abexp, aes(y = model_name, x = auc, fill = model_name)) +
  geom_boxplot() +
  scale_fill_manual(values = colors_here) +
  stat_compare_means(
    comparisons = list(c("AbExp \n(AbSplice2)", "AbExp \n(AbSplice1)")),
    method = "wilcox.test",
    paired = TRUE,
    method.args = list(alternative = "greater", exact = FALSE)
  ) +
  labs(y = NULL, x = "auPRC", title = "") +
  # theme_cowplot(font_size = fontsize) +
  theme_cowplot() +
  theme(legend.position = "none")

# Combine with ggarrange
abexp_performance2 <- ggarrange(
  plot_absplice,
  plot_abexp,
  ncol = 1,
  align = "v",
  labels = NULL
)
abexp_performance2

title_plot <- ggdraw() +
  draw_label(
    label = "Aberrant underexpression prediction\n(GTEx, across tissues)",
    x = 0,               # far left
    hjust = 0,           # align to left edge
    fontface = "bold",
    size = fontsize*1.142857
  )
abexp_performance3 <- ggarrange(
  title_plot,
  plot_absplice,
  plot_abexp,
  ncol = 1,
  heights = c(0.1, 1, 1),
  align = "v"
)
abexp_performance3


# ========== Assemble Fig. 2
fig_2 <- ggarrange(
  ggarrange(
    fig_pr_curve_gtex, 
    fig_box, 
    ncol = 2, labels = c("a", "b")),
  ggarrange(
    fig_pr_curve_gvex, 
    abexp_performance3,
    ncol = 2, labels = c("c", "d")),
  nrow=2, heights = c(1, 1)
)
fig_2

ggsave(
  filename = file.path(DATA_DIR, "AbSplice2_figures/Fig_2.svg"),
  plot = fig_2,
  width = 30,
  height = 25,
  units = "cm",
  device = "svg"
)





