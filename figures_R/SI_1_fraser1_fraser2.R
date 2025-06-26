source("config.R")

# ===================== FRASER1 -> FRASER2 =====================

df <- fread(pr_curve_FRASER1_and_FRASER2_ground_truth)
color_list = absplice_model_colors[unique(df$model)]

# Specify the order you want for the legend:
desired_order <- c(
  "AbSplice-DNA (trained on FRASER1, evaluated on FRASER1)",
  'AbSplice-DNA (trained on FRASER2, evaluated on FRASER1)',
  'AbSplice-DNA (trained on FRASER1, evaluated on FRASER2)',
  'AbSplice-DNA (trained on FRASER2, evaluated on FRASER2)'
)
# Convert 'model' to factor with the specified levels:
df$model <- factor(df$model, levels = rev(desired_order))


plot_pr_curve <- function(df, color_list, ylim=0.2, breaks_x=0.2, minor_breaks_x=0.2, breaks_y=0.1, minor_breaks_y=0.02, title='All tissues') {
  g <- (
    ggplot(df, aes(x=recall, y=precision, color=model))
    + geom_step(direction = "vh", linewidth=1)
    + scale_color_manual(values=color_list, guide=guide_legend(reverse = TRUE, nrow=4))
    + theme_cowplot(font_size = fontsize)
    + scale_x_continuous(limits = c(-0.01, 1.05), breaks = seq(0, 1, breaks_x), minor_breaks = seq(0, 1, minor_breaks_x))
    + scale_y_continuous(limits = c(-0.01, ylim), breaks = seq(0, 1, breaks_y), minor_breaks = seq(0, 1, minor_breaks_y))
    + labs(x='Recall', y='Precision', color='', title=title)
    + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    + theme(
      legend.position = "none",
    )
  )
  return(g)
}

fig <- plot_pr_curve(
  df,
  ylim=0.4,
  color_list,
  title='All tissues',
)
fig

df <- fread(pr_curve_FRASER1_and_FRASER2_ground_truth_boxplot)
df <- process_models_box_plot(df, unique(df$model))

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
    ggplot(df, aes(x=`Average Precision Score`, y=model_label, fill=model_label))
    + scale_fill_manual(values=color_list, guide=guide_legend(reverse = FALSE, nrow=4))
    + geom_boxplot()
    + xlim(0, x_lim)
    + theme_cowplot(font_size = fontsize)
    + theme(legend.title = element_blank())
    + theme(legend.position = "bottom")
    + labs(
      x="auPRC",
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
    + theme(
      axis.title.x = element_blank(),
      axis.text.x  = element_blank(),
    )
  )
  return(g)
}

label_dict <- c(
  "AbSplice-DNA (trained on FRASER1, evaluated on FRASER1)" = "AbSplice1 (trained on FRASER1, evaluated on FRASER1)",
  "AbSplice-DNA (trained on FRASER2, evaluated on FRASER1)" = "AbSplice1 (trained on FRASER2, evaluated on FRASER1)",
  "AbSplice-DNA (trained on FRASER1, evaluated on FRASER2)" = "AbSplice1 (trained on FRASER1, evaluated on FRASER2)",
  "AbSplice-DNA (trained on FRASER2, evaluated on FRASER2)" = "AbSplice1 (trained on FRASER2, evaluated on FRASER2)"
)

df$model_label <- factor(label_dict[as.character(df$model)], levels = label_dict[desired_order])
color_list <- setNames(absplice_model_colors[desired_order], label_dict[desired_order])

fig_box <- plot_box_plot(
  df,
  color_list, 
  coord_flip=TRUE,
  title='Across tissues'
)
fig_box


df <- fread(feature_importance_absplice1_fraser1_fraser2)

color_dict <- c(
  "fraser1"= "red",
  "fraser2"= "orange"
)

facet_labels <- c(
  "delta_logit_psi"= "MMSplice SpliceMap (ΔLogit Ψ)",
  "delta_score"= "SpliceAI",
  "delta_psi"= "MMSplice SpliceMap (ΔΨ)",
  "splice_site_is_expressed"= "Splice Site Expressed (median read cutoff)"
)

x_labels <- c(
  "delta_logit_psi" = "ΔLogit Ψ",
  "delta_score" = "Delta Score",
  "delta_psi" = "ΔΨ",
  "splice_site_is_expressed" = "Median split read counts of splice site"
)

# Define font sizes as variables
fontsize_axis_labels <- 14
fontsize_facet_labels <- 14
fontsize_legend <- 14

fig_features <- ggplot(df, aes(x = score_value, y = `Contribution Score`, color = score)) +
  geom_step(direction='vh') +
  facet_wrap(~feature_name, ncol = 1, scales = "free_x",
             labeller = labeller(feature_name = facet_labels)) +
  scale_color_manual(
    values = color_dict,
    labels = c("fraser1" = "FRASER1", "fraser2" = "FRASER2"),
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

fig_2 <- ggarrange(
  ggarrange(fig, fig_box, ncol = 1, heights = c(5,5), labels = c("a", "b")),  # Left column: fig and fig_box stacked
  fig_features,  # Right column: fig_features
  ncol = 2, widths = c(2,2), labels = c("", "c")  # Arrange columns, adjust widths as needed
)
fig_2

ggsave(
  filename = file.path(DATA_DIR, "AbSplice2_figures/SI_1_fraser1_fraser2.svg"),
  plot = fig_2,
  width = 30,
  height = 25,
  units = "cm",
  device = "svg"
)
