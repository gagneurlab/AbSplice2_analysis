source("config.R")

# ===================== Singed features =====================

plot_pr_curve <- function(df, color_list, ylim=0.2, breaks_x=0.2, minor_breaks_x=0.2, breaks_y=0.1, minor_breaks_y=0.02, title='All tissues') {
  g <- (
    ggplot(df, aes(x=recall, y=precision, color=model_label))
    + geom_step(direction = "vh", linewidth=1.2)
    + scale_color_manual(values=color_list, guide=guide_legend(reverse = FALSE, nrow=4))
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

colors_here <- c(
  # Pangolin-Based Models (Expanded Orange Range)
  "Pangolin Score" = "#ffd699",   # Light Orange (Base)
  "Pangolin max score" = "#ffd699",   # Light Orange (Base)
  "AbSplice-DNA 1 (Pangolin max score)" = "#ffac5e",  # Medium Light Orange
  "AbSplice-DNA 1 (Pangolin max score, signed)" = "#ff7f0e", # Medium Orange
  "AbSplice-DNA 1 (Pangolin all scores, signed)" = "#a84300", # Deepest Orange
  
  # SpliceAI-Based Models (Expanded Blue Range)
  "SpliceAI Delta Score" = "#b3d9ff", # Light Blue (Base)
  "SpliceAI max score" = "#b3d9ff", # Light Blue (Base)
  "AbSplice-DNA 1 (SpliceAI max score)" = "#72b0e0",  # Medium Light Blue
  "AbSplice-DNA 1 (SpliceAI max score, signed)"= "#1f77b4", # Medium Blue
  "AbSplice-DNA 1 (SpliceAI all scores, signed)" = "#08306b" # Deepest Blue
)

df <- fread(pr_curve_signed_features_SpliceAI_Pangolin)

chosen_models <- c(
  "AbSplice-DNA 1 (Pangolin max score)",
  "AbSplice-DNA 1 (Pangolin max score, signed)",
  "AbSplice-DNA 1 (Pangolin all scores, signed)"
)
df_here <- df[model %in% chosen_models]
df_here <- process_models_pr_curve(df_here, chosen_models)

color_list = colors_here[unique(df_here$model)]

label_dict <- c(
  "AbSplice-DNA 1 (Pangolin max score)" = "AbSplice1 (Pangolin max score)",
  "AbSplice-DNA 1 (Pangolin max score, signed)" = "AbSplice1 (Pangolin max score, signed)",
  "AbSplice-DNA 1 (Pangolin all scores, signed)" = "AbSplice1 (Pangolin all scores, signed)"
)

df_here$model_label <- factor(label_dict[as.character(df_here$model)], levels = label_dict[chosen_models])
color_list <- setNames(colors_here[names(label_dict)], label_dict[names(label_dict)])

fig1 <- plot_pr_curve(
  df_here,
  ylim=0.5,
  color_list,
  title='',
)
fig1

chosen_models <- c(
  "AbSplice-DNA 1 (SpliceAI max score)",
  "AbSplice-DNA 1 (SpliceAI max score, signed)",
  "AbSplice-DNA 1 (SpliceAI all scores, signed)"
)
df_here <- df[model %in% chosen_models]
df_here <- process_models_pr_curve(df_here, chosen_models)

color_list = colors_here[unique(df_here$model)]

label_dict <- c(
  "AbSplice-DNA 1 (SpliceAI max score)" = "AbSplice1 (SpliceAI max score) := AbSplice1",
  "AbSplice-DNA 1 (SpliceAI max score, signed)" = "AbSplice1 (SpliceAI max score, signed)",
  "AbSplice-DNA 1 (SpliceAI all scores, signed)" = "AbSplice1 (SpliceAI all scores, signed)"
)

df_here$model_label <- factor(label_dict[as.character(df_here$model)], levels = label_dict[chosen_models])
color_list <- setNames(colors_here[names(label_dict)], label_dict[names(label_dict)])


fig2 <- plot_pr_curve(
  df_here,
  ylim=0.5,
  color_list,
  title='',
)
fig2


chosen_models <- c(
  "SpliceAI max score",
  "Pangolin max score",
  "AbSplice-DNA 1 (SpliceAI all scores, signed)",
  "AbSplice-DNA 1 (Pangolin all scores, signed)"
)

df <- df %>%
  mutate(model = case_when(
    model == "SpliceAI Delta Score" ~ "SpliceAI max score",
    model == "Pangolin Score" ~ "Pangolin max score",
    TRUE ~ model  # keep other values unchanged
  ))

df_here <- df[model %in% chosen_models]
df_here <- process_models_pr_curve(df_here, chosen_models)

color_list = colors_here[unique(df_here$model)]

label_dict <- c(
  "SpliceAI max score" = "SpliceAI max score",
  "Pangolin max score" = "Pangolin max score",
  "AbSplice-DNA 1 (SpliceAI all scores, signed)" = "AbSplice1 (SpliceAI all scores, signed)",
  "AbSplice-DNA 1 (Pangolin all scores, signed)" = "AbSplice1 (Pangolin all scores, signed)"
)

df_here$model_label <- factor(label_dict[as.character(df_here$model)], levels = label_dict[chosen_models])
color_list <- setNames(colors_here[names(label_dict)], label_dict[names(label_dict)])

fig3 <- plot_pr_curve(
  df_here,
  ylim=0.5,
  color_list,
  title='',
)
fig3

fig_all <- ggarrange(
  ncol = 3, widths = c(2,2), heights = c(5,5), labels = c("a", "b", "c"),
  fig2, fig1, fig3
)
fig_all

ggsave(
  filename = file.path(DATA_DIR, "AbSplice2_figures/SI_2_signed_features.svg"),
  plot = fig_all,
  width = 40,
  height = 15,
  units = "cm",
  device = "svg"
)

