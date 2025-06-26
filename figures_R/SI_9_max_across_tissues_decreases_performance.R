source("config.R")

# ===================== Max across tissues =====================

plot_pr_curve <- function(df, color_list, ylim=0.2, breaks_x=0.2, minor_breaks_x=0.2, breaks_y=0.1, minor_breaks_y=0.02, title='All tissues') {
  g <- (
    ggplot(df, aes(x=recall, y=precision, color=model_label))
    + geom_step(direction = "vh", linewidth=0.9)
    + scale_color_manual(values=color_list, guide=guide_legend(reverse = TRUE, nrow=3))
    + theme_cowplot(font_size = fontsize)
    + scale_x_continuous(limits = c(-0.01, 1.05), breaks = seq(0, 1, breaks_x), minor_breaks = seq(0, 1, minor_breaks_x))
    + scale_y_continuous(limits = c(-0.01, ylim), breaks = seq(0, 1, breaks_y), minor_breaks = seq(0, 1, minor_breaks_y))
    + labs(x='Recall', y='Precision', color='', title=title)
    + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    + theme(
      legend.position = "bottom",
    )
  )
  return(g)
}

df <- fread(pr_curve_max_aggregation_decreases_performance)

colors_here <- c(
  # Base Models
  "SpliceAI" = "#72b0e0",  # Medium Light Blue
  "Pangolin" = "#ffac5e",  # Medium Light Orange
  
  # Core AbSplice-DNA Models (Darkest)
  "AbSplice-DNA 1" = "#67000d",  # Darkest Red
  "AbSplice-DNA 2" = "#67000d",  # Darkest Red (Same as DNA 1)
  
  # Max Main Tissues (Slightly Lighter)
  "AbSplice DNA 1 (max main tissues)" = "#a50f15",  # Deep Red
  "AbSplice DNA 2 (max main tissues)" = "#a50f15",  # Deep Red (Same as DNA 1)
  
  # Max All (Lightest)
  "AbSplice DNA 1 (max all)" = "#ef3b2c",  # Lightest Red
  "AbSplice DNA 2 (max all)" = "#ef3b2c"   # Lightest Red (Same as DNA 1)
)

chosen_models <- c(
  "SpliceAI",
  "Pangolin",
  "AbSplice DNA 1 (max all)",
  "AbSplice DNA 1 (max main tissues)",
  "AbSplice-DNA 1"
)
df_here <- df[model %in% chosen_models]
df_here <- df_here[outlier %in% "FRASER1"]
df_here <- process_models_pr_curve(df_here, chosen_models)

color_list = colors_here[unique(df_here$model)]

label_dict <- c(
  "SpliceAI" = "SpliceAI",
  "Pangolin" = "Pangolin",
  "AbSplice DNA 1 (max all)" = "AbSplice1 (max over all tissues)",
  "AbSplice DNA 1 (max main tissues)" = "AbSplice1 (max over main tissues)",
  "AbSplice-DNA 1" = "AbSplice1"
)

df_here$model_label <- factor(label_dict[as.character(df_here$model)], levels = label_dict[chosen_models])
color_list <- setNames(colors_here[names(label_dict)], label_dict[names(label_dict)])

fig1 <- plot_pr_curve(
  df_here,
  ylim=1.0,
  color_list,
  title='Evaluated on FRASER1 outliers',
)
fig1


chosen_models <- c(
  "SpliceAI",
  "Pangolin",
  "AbSplice DNA 2 (max all)",
  "AbSplice DNA 2 (max main tissues)",
  "AbSplice-DNA 2"
)
df_here <- df[model %in% chosen_models]
df_here <- df_here[outlier %in% "FRASER2"]
df_here <- process_models_pr_curve(df_here, chosen_models)

color_list = colors_here[unique(df_here$model)]

label_dict <- c(
  "SpliceAI" = "SpliceAI",
  "Pangolin" = "Pangolin",
  "AbSplice DNA 2 (max all)" = "AbSplice2 (max over all tissues)",
  "AbSplice DNA 2 (max main tissues)" = "AbSplice2 (max over main tissues)",
  "AbSplice-DNA 2" = "AbSplice2"
)

df_here$model_label <- factor(label_dict[as.character(df_here$model)], levels = label_dict[chosen_models])
color_list <- setNames(colors_here[names(label_dict)], label_dict[names(label_dict)])

fig2 <- plot_pr_curve(
  df_here,
  ylim=1.0,
  color_list,
  title='Evaluated on FRASER2 outliers',
)
fig2


chosen_models <- c(
  "SpliceAI",
  "Pangolin",
  "AbSplice DNA 2 (max all)",
  "AbSplice DNA 2 (max main tissues)",
  "AbSplice-DNA 2"
)
df_here <- df[model %in% chosen_models]
df_here <- df_here[outlier %in% "FRASER1"]
df_here <- process_models_pr_curve(df_here, chosen_models)

color_list = colors_here[unique(df_here$model)]

label_dict <- c(
  "SpliceAI" = "SpliceAI",
  "Pangolin" = "Pangolin",
  "AbSplice DNA 2 (max all)" = "AbSplice2 (max over all tissues)",
  "AbSplice DNA 2 (max main tissues)" = "AbSplice2 (max over main tissues)",
  "AbSplice-DNA 2" = "AbSplice2"
)

df_here$model_label <- factor(label_dict[as.character(df_here$model)], levels = label_dict[chosen_models])
color_list <- setNames(colors_here[names(label_dict)], label_dict[names(label_dict)])

fig3 <- plot_pr_curve(
  df_here,
  ylim=1.0,
  color_list,
  title='Evaluated on FRASER1 outliers',
)
fig3

chosen_models <- c(
  "SpliceAI",
  "Pangolin",
  "AbSplice DNA 1 (max all)",
  "AbSplice DNA 1 (max main tissues)",
  "AbSplice-DNA 1"
)
df_here <- df[model %in% chosen_models]
df_here <- df_here[outlier %in% "FRASER2"]
df_here <- process_models_pr_curve(df_here, chosen_models)

color_list = colors_here[unique(df_here$model)]

label_dict <- c(
  "SpliceAI" = "SpliceAI",
  "Pangolin" = "Pangolin",
  "AbSplice DNA 1 (max all)" = "AbSplice1 (max over all tissues)",
  "AbSplice DNA 1 (max main tissues)" = "AbSplice1 (max over main tissues)",
  "AbSplice-DNA 1" = "AbSplice1"
)

df_here$model_label <- factor(label_dict[as.character(df_here$model)], levels = label_dict[chosen_models])
color_list <- setNames(colors_here[names(label_dict)], label_dict[names(label_dict)])

fig4 <- plot_pr_curve(
  df_here,
  ylim=1.0,
  color_list,
  title='Evaluated on FRASER2 outliers',
)
fig4

fig_all <- ggarrange(
  fig1, fig4, fig3, fig2, labels=c('a', 'b', 'c','d'), # Arrange all four figures
  nrow = 2, ncol = 2       # Set rows and columns
)
fig_all

ggsave(
  filename = file.path(DATA_DIR, "AbSplice2_figures/SI_5_max_across_tissues_decreases_performance.svg"),
  plot = fig_all,
  width = 30,
  height = 25,
  units = "cm",
  device = "svg"
)


