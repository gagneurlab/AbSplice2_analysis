source("config.R")

plot_pr_curve <- function(df, color_list, ylim=0.2, breaks_x=0.2, minor_breaks_x=0.2, breaks_y=0.1, minor_breaks_y=0.02, title='All tissues') {
  g <- (
    ggplot(df, aes(x=recall, y=precision, color=model))
    + geom_step(direction = "vh", linewidth=0.9)
    + scale_color_manual(values = color_list, 
                         labels = c("AbSplice-DNA (SpliceAI -> Pangolin)" = "AbSplice-DNA with binary splice site usage", 
                                    "AbSplice-DNA 2" = "AbSplice-DNA with continuous splice site usage"),
                         guide = guide_legend(reverse = TRUE, nrow = 4))
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

df <- fread(pr_curve_final_absplice2_model)

colors_here <- c(
  "AbSplice-DNA (SpliceAI -> Pangolin)" = "#2ca02c",  # Green
  "AbSplice-DNA 2" = "#d62728"  # Strong Red
)
chosen_models <- c(
  "AbSplice-DNA (SpliceAI -> Pangolin)",
  "AbSplice-DNA 2"
)

df_here <- df[model %in% chosen_models]
df_here <- process_models_pr_curve(df_here, chosen_models)
df_here$model <- factor(df_here$model, levels = chosen_models)

color_list = colors_here[unique(df_here$model)]

pr_curve <- plot_pr_curve(
  df_here,
  ylim=1.0,
  color_list,
  title='Aberrant splicing prediction \n(GTEx, all tissues)',
)
pr_curve


## ========== feature importance
df <- fread(file.path(DATA_DIR, 'feature_importance_binary_and_continous_splice_site_usage.csv'))

color_dict <- c(
  "continous"= "red",
  "binary"= "orange"
)

facet_labels <- c(
  "delta_logit_psi"= "MMSplice SpliceMap (ΔLogit Ψ)",
  "delta_score"= "SpliceAI",
  "delta_psi"= "MMSplice SpliceMap (ΔΨ)",
  "median_n"= "Splice site usage (median read counts)"
)

x_labels <- c(
  "delta_logit_psi" = "ΔLogit Ψ",
  "delta_score" = "Delta Score",
  "delta_psi" = "ΔΨ",
  "median_n" = "Median split read counts of splice site"
)

# Define font sizes as variables
fontsize_axis_labels <- 14
fontsize_facet_labels <- 14
fontsize_legend <- 14

fig_features_all <- ggplot(df, aes(x = score_value, y = `Contribution Score`, color = model)) +
  geom_step(direction='vh') +
  facet_wrap(~feature, ncol = 1, scales = "free_x",
             labeller = labeller(feature = facet_labels)) +
  scale_color_manual(
    values = color_dict,
    labels = c("continous" = "continuous", "binary" = "binary"),
    guide=guide_legend(reverse = FALSE, nrow=2)
  ) +
  theme_minimal() +
  labs(
    x = "Score Value",
    y = "Contribution Score",
    color = "Model trained on"
  ) +
  theme_cowplot(font_size = fontsize) +
  theme(
    plot.margin = margin(5, 5, 5, 5),
    text = element_text(size = fontsize_axis_labels), # Default text size (can be adjusted further)
    strip.text = element_text(size = fontsize_facet_labels),  # Facet labels fontsize
    axis.title.x = element_text(size = fontsize_axis_labels),  # X-axis label fontsize
    axis.title.y = element_text(size = fontsize_axis_labels),  # Y-axis label fontsize
    axis.text.x = element_text(size = 10), # Keep axis text smaller for readability
    legend.position = "bottom",               # Move legend to the top
    legend.text = element_text(size = fontsize_legend),   # Legend text fontsize
    legend.title = element_text(size = fontsize_legend),  # Legend title fontsize
    strip.background = element_rect(colour="white", fill="white")
  )

fig_features_all

fig_features_median_n <- ggplot(df[feature %in% c('median_n')], aes(x = score_value, y = `Contribution Score`, color = model)) +
  geom_step(direction='vh') +
  scale_color_manual(
    values = color_dict,
    labels = c("continous" = "continuous splice site usage", "binary" = "binary splice site usage"),
    guide=guide_legend(reverse = FALSE, nrow=2)
  ) +
  theme_minimal() +
  labs(
    x = "Median split read counts of splice site",
    y = "Contribution Score",
    color = "Model trained on"
  ) +
  theme_cowplot(font_size = fontsize) +
  theme(
    plot.margin = margin(5, 5, 5, 5),
    text = element_text(size = fontsize_axis_labels), # Default text size (can be adjusted further)
    strip.text = element_text(size = fontsize_facet_labels),  # Facet labels fontsize
    axis.title.x = element_text(size = fontsize_axis_labels),  # X-axis label fontsize
    axis.title.y = element_text(size = fontsize_axis_labels),  # Y-axis label fontsize
    axis.text.x = element_text(size = 10), # Keep axis text smaller for readability
    legend.position = "bottom",               # Move legend to the top
    legend.text = element_text(size = fontsize_legend),   # Legend text fontsize
    legend.title = element_text(size = fontsize_legend),  # Legend title fontsize
    strip.background = element_rect(colour="white", fill="white")
  ) +
  theme(
    legend.position = c(0.95, 0.3),
    legend.justification = c("right", "bottom")
  )

fig_features_median_n


## ============= absplice median_n vs binary

df_plot <- fread(file.path(DATA_DIR, 'median_n_reranking_absplice_scores.csv'))

# Define axis limits
min_val <- min(df_plot$AbSplice_binary, df_plot$AbSplice_continous)
max_val <- max(df_plot$AbSplice_binary, df_plot$AbSplice_continous)
max_val <- 0.6  # Adjusted manually

# Set max count value for color mapping
max_count <- 1000  # Values above this will be shown in the same color

df_plot$median_n_category <- factor(df_plot$median_n_category, levels = c("<10", "10-100", ">100"))

# Define label dictionaries
facet_labels_ytest <- c("0" = "non-outliers", "1" = "outliers")
facet_labels_median <- c(
  "<10" = "<10",
  "10-100" = "10–100",
  ">100" = ">100"
)
abplice_binary_vs_continuous <- ggplot(df_plot, aes(
  x = AbSplice_binary, 
  y = AbSplice_continous
)) +
  geom_bin2d(bins = 50) +
  scale_fill_viridis_c(
    trans = "log10",
    limits = c(1, max_count),
    breaks = c(1, 10, 100, 1000),
    labels = c("1", "10", "100", "1000+"),
    option = "viridis"
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  facet_grid(
    y_test ~ median_n_category, 
    labeller = labeller(
      y_test = facet_labels_ytest,
      median_n_category = facet_labels_median
    )
  ) +
  coord_cartesian(xlim = c(min_val, max_val), ylim = c(min_val, max_val)) +
  labs(
    x = "AbSplice-DNA with binary splice site usage", 
    y = "AbSplice-DNA with continuous splice site usage"
  ) +
  theme_minimal() +
  theme(
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9),
    axis.title = element_text(size = 10)
  ) +
  theme_cowplot(font_size = fontsize)

abplice_binary_vs_continuous


fig <- ggarrange(
  ggarrange(pr_curve, fig_features_median_n, nrow = 2, labels = c("a", "b")), 
  ggarrange(abplice_binary_vs_continuous, ncol = 1, heights = c(5,5), labels = c("c")),
  ncol=2, widths=c(2,3))
fig

ggsave(
  filename = file.path(DATA_DIR, "AbSplice2_figures/SI_3_continuous_median_n.svg"),
  plot = fig,
  width = 35,
  height = 25,
  units = "cm",
  device = "svg"
)

