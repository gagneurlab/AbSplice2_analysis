setwd("/data/nasif12/home_if12/wagnern/Projects/gitlab_gagneurlab/AbSplice_analysis/workflow/scripts/figures_R")
source("config.R")

# ======== Maximum across all junctions
# Load data
df_all_save <- fread('/s/project/absplice/AbSplice_analysis/data/results/source_data/GEL_example_df_all_chr8_38285617_C_T.csv')

# Select tissue
tissue_select <- "brain"
df_plot <- df_all_save[df_all_save$tissue_type == tissue_select, ]

# Set desired timepoint order
tp_description_order <- c(
  "4 wpc", "5 wpc", "6 wpc", "7 wpc", "8 wpc",
  "9–11 wpc", "12 wpc", "13 wpc", "16 wpc",
  "18–19 wpc", "20 wpc", "Newborn–Toddler",
  "School–Young Adult", "Young Mid Age", "Older Mid–Senior", "GTEx"
)

# Convert tp_description to factor with correct order
df_plot$tp_description <- factor(df_plot$tp_description, levels = tp_description_order)

df_plot <- df_plot[df_plot$tp_description != "GTEx", ]


# Plot
GEL_case_all <- ggplot(df_plot, aes(x = tp_description, y = model_27, color = model_27_category)) +
  geom_point() +
  scale_color_manual(values = c("high" = "red", "medium" = "orange", "low" = "blue")) +
  labs(
    y = "AbSplice2",
    x = "",
    color = "Score Category",
    title = 'Maximum variant effect'
  ) +
  scale_y_continuous(limits = c(0, max(df_plot$model_27, na.rm = TRUE))) +
  theme_cowplot(font_size = fontsize) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  theme(
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top")
  )

GEL_case_all


# ========= affected junction

# Load data
df_save <- fread('/s/project/absplice/AbSplice_analysis/data/results/source_data/GEL_example_df_chr8_38285617_C_T.csv')  # update the path if needed

# Set selections
tissue_select <- "brain"
junction_select <- "chr8:38285615-38285863:-"

# Filter and prepare
df_plot <- df_save[df_save$tissue_type == tissue_select & df_save$junction == junction_select, ]
df_plot$tp_description <- factor(df_plot$tp_description, levels = tp_description_order)

# Filter out GTEx
df_plot <- df_plot[df_plot$tp_description != "GTEx", ]

# Plot
GEL_weak_site <- ggplot(df_plot, aes(x = tp_description, y = model_27, color = model_27_category)) +
  geom_point() +
  scale_color_manual(values = c("high" = "red", "medium" = "orange", "low" = "blue")) +
  labs(
    y = "AbSplice2",
    x = "",
    # title = paste("Variant effect on weak splice site:", junction_select),
    title = paste("Variant effect on weak splice site"),
    color = "Score Category"
  ) +
  scale_y_continuous(limits = c(0, max(df_plot$model_27, na.rm = TRUE))) +
  theme_cowplot(font_size = fontsize) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  theme(
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top")
  )
GEL_weak_site


placeholder <- ggplot() + 
  theme_void() +
  theme(panel.background = element_rect(fill = NA, color = NA)) 

fig_GEL <- ggarrange(
  placeholder, labels=c('a'),
  ggarrange(
    GEL_case_all,
    GEL_weak_site,
    nrow=2, labels=c('b', 'c')
  ), ncol=2, widths=c(2,1)
)
fig_GEL


ggsave(
  filename = "/s/project/absplice/AbSplice_analysis/data/results/source_data/AbSplice2_figures/Fig_4_AbSplice_preds.svg",
  plot = fig_GEL,
  width = 35,
  height = 10,
  units = "cm",
  device = "svg"
)






# ================
# Load data
df_all_save <- fread('/s/project/absplice/AbSplice_analysis/data/results/source_data/GEL_example_df_all_chr8_38285617_C_T.csv')
df_save <- fread('/s/project/absplice/AbSplice_analysis/data/results/source_data/GEL_example_df_chr8_38285617_C_T.csv')

# Settings
tissue_select <- "brain"
junction_select <- "chr8:38285615-38285863:-"
tp_description_order <- c(
  "4 wpc", "5 wpc", "6 wpc", "7 wpc", "8 wpc",
  "9–11 wpc", "12 wpc", "13 wpc", "16 wpc",
  "18–19 wpc", "20 wpc", "Newborn–Toddler",
  "School–Young Adult", "Young Mid Age", "Older Mid–Senior", "GTEx"
)

# Prepare df1: maximum across all junctions
df_all <- df_all_save %>%
  filter(tissue_type == tissue_select, tp_description != "GTEx") %>%
  mutate(
    tp_description = factor(tp_description, levels = tp_description_order),
    facet_label = "Maximum variant effect"
  )

# Prepare df2: affected junction only
df_junction <- df_save %>%
  filter(tissue_type == tissue_select, junction == junction_select, tp_description != "GTEx") %>%
  mutate(
    tp_description = factor(tp_description, levels = tp_description_order),
    facet_label = "Variant effect on weak splice site"
  )

# Combine both
df_combined <- bind_rows(df_all, df_junction)

# Plot
combined_plot <- ggplot(df_combined, aes(x = tp_description, y = model_27, color = model_27_category)) +
  geom_point() +
  scale_color_manual(values = c("high" = "red", "medium" = "orange", "low" = "blue")) +
  facet_wrap(~facet_label, ncol = 1, scales = "fixed") +
  # geom_hline(yintercept = 0.1, linetype = "dashed", color = "blue") +
  # geom_hline(yintercept = 0.2, linetype = "dashed", color = "orange") +
  labs(
    x = "",
    y = "AbSplice2",
    color = "Score Category"
  ) +
  scale_y_continuous(limits = c(0, max(df_combined$model_27, na.rm = TRUE))) +
  theme_cowplot(font_size = fontsize) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = fontsize)
  )  
  # theme(
  #   legend.position = c(0.95, 0.95),
  #   legend.justification = c("right", "top")
  # )
combined_plot

fig_GEL2 <- ggarrange(
  placeholder, combined_plot, labels=c('a', 'b'), ncol=2, widths=c(2,1)
)
fig_GEL2


ggsave(
  filename = "/s/project/absplice/AbSplice_analysis/data/results/source_data/AbSplice2_figures/Fig_4_AbSplice_preds.svg",
  plot = combined_plot,
  # width = 12,
  # height = 12,
  width = 20,
  height = 15,
  units = "cm",
  device = "svg"
)

