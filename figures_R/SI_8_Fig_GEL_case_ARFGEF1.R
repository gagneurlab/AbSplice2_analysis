source("config.R")

# Load data
df_all_save <- fread(file.path(DATA_DIR, 'GEL_example_df_all_chr8_68115487_C_T.csv'))
df_save <- fread(file.path(DATA_DIR, 'GEL_example_df_chr8_68115487_C_T.csv'))

# Settings
tissue_select <- "brain"
junction_select <- "chr8:68115484-68116914:-"


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
combined_plot_ARFGEF1 <- ggplot(df_combined, aes(x = tp_description, y = model_27, color = model_27_category)) +
  geom_point() +
  scale_color_manual(values = c("high" = "red", "medium" = "orange", "low" = "blue")) +
  facet_wrap(~facet_label, ncol = 1, scales = "fixed") +
  labs(
    title = "ARFGEF1",
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
combined_plot_ARFGEF1


# FN1
df_all <- fread(file.path(DATA_DIR, 'GEL_example_df_all_chr2_216237105_G_C.csv'))

# Timepoint order
tp_description_order <- c(
  "4 wpc", "5 wpc", "6 wpc", "7 wpc", "8 wpc",
  "9–11 wpc", "12 wpc", "13 wpc", "16 wpc",
  "18–19 wpc", "20 wpc", "Newborn–Toddler",
  "School–Young Adult", "Young Mid Age", "Older Mid–Senior", "GTEx"
)

# Preprocess
df_all <- df_all %>%
  filter(tp_description %in% tp_description_order & tp_description != "GTEx") %>%
  mutate(
    tp_description = factor(tp_description, levels = tp_description_order),
    tissue_type = factor(tissue_type)
  )

# Plot
plot_facet_FN1 <- ggplot(df_all, aes(x = tp_description, y = model_27)) +
  geom_point(size = 1) +
  facet_wrap(~ tissue_type) +
  labs(
    title = "FN1",
    y = "AbSplice2",
    x = ""
  ) +
  scale_y_continuous(limits = c(0, max(df_all$model_27, na.rm = TRUE))) +
  theme_cowplot(font_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

plot_facet_FN1

GEL_cases <- ggarrange(
  combined_plot_ARFGEF1,
  plot_facet_FN1,
  ncol=2, widths=c(1,1.5), labels=c('a', 'b')
)
GEL_cases


ggsave(
  filename = file.path(DATA_DIR, "AbSplice2_figures/SI_Fig_GEL_cases_AbSplice_preds.svg"),
  plot = GEL_cases,
  width = 40,
  height = 25,
  units = "cm",
  device = "svg"
)


