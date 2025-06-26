source("config.R")

library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

# =========== dev example region
dev_example <- readPNG(file.path(DATA_DIR, "fig3_only_example.png"))

dev_example_fig <- as.raster(dev_example)
dev_example_fig <- rasterGrob(dev_example_fig, interpolate = FALSE)
dev_example_fig

# # ================================================  Developmental variant
df <- fread(file.path(DATA_DIR, 'developmental_variants.csv'))

# Define timepoint order and relabeling
timepoints <- c(paste0("t", 1:15), "t_gtex")
timepoint_labels <- c(
  t1 = "4 wpc",
  t2 = "5 wpc",
  t3 = "6 wpc",
  t4 = "7 wpc",
  t5 = "8 wpc",
  t6 = "9–11 wpc",
  t7 = "12 wpc",
  t8 = "13 wpc",
  t9 = "16 wpc",
  t10 = "18–19 wpc",
  t11 = "20 wpc",
  t12 = "Newborn–Toddler",
  t13 = "School–Young Adult",
  t14 = "Young Mid Age",
  t15 = "Older Mid–Senior",
  t_gtex = "GTEx"
)

df$timepoint <- factor(df$timepoint, levels = timepoints, labels = timepoint_labels)

# Add AbSplice2 category for coloring points
df <- df %>%
  mutate(absplice_category = case_when(
    AbSplice2 > 0.2 ~ "high",
    AbSplice2 > 0.1 ~ "medium",
    AbSplice2 > 0.05 ~ "low",
    TRUE ~ "none"
  ))

# Reshape to long format
df_long <- df %>%
  pivot_longer(cols = c(ref_psi, AbSplice2),
               names_to = "metric",
               values_to = "value") %>%
  mutate(
    absplice_category = ifelse(metric == "AbSplice2", absplice_category, NA),
    # Create facet labels
    facet_metric = recode(metric,
                          ref_psi = "Reference PSI",
                          AbSplice2 = "AbSplice2"),
    facet_metric = factor(facet_metric, levels = c("Reference PSI", "AbSplice2"))
  )

# Define color mapping
category_colors <- c(
  "high" = "red",
  "medium" = "orange",
  "low" = "blue",
  "none" = "grey70"
)

# Horizontal cutoff lines only for AbSplice2 facet
cutoffs_df <- data.frame(
  facet_metric = factor("AbSplice2", levels = c("Reference PSI", "AbSplice2")),
  # yintercept = c(0.01, 0.05, 0.2),
  yintercept = c(0.05, 0.1, 0.2),
  cutoff_level = factor(c("low", "medium", "high"), levels = c("high", "medium", "low")),
  color = c("blue", "orange", "red")
)

df_long <- df_long[df_long$timepoint != "GTEx", ]


# Final plot
dev_var_plot <- ggplot(df_long, aes(x = timepoint, y = value)) +
  # AbSplice2 colored points
  geom_point(
    data = df_long %>% filter(metric == "AbSplice2"),
    aes(color = absplice_category),
    size = 2,
    show.legend = FALSE
  ) +
  # Horizontal cutoff lines for AbSplice2
  geom_hline(
    data = cutoffs_df,
    aes(yintercept = yintercept, color = cutoff_level),
    linetype = "dotted",
    show.legend = FALSE
  ) +
  # ref_psi black points
  geom_point(
    data = df_long %>% filter(metric == "ref_psi"),
    color = "black",
    size = 2
  ) +
  scale_color_manual(
    values = category_colors,
    name = "Cutoffs",
    breaks = c("high", "medium", "low")
  ) +
  theme_minimal(base_size = fontsize) +
  theme_cowplot(font_size = fontsize) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    # title = "Dynamic Metrics",
    title = "",
    y = "Value",
    x = "Timepoint"
  ) +
  facet_wrap(vars(facet_metric), scales = "free_y", ncol = 1)

# Show plot
dev_var_plot


# ================================================  Static scores
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

# Step 1: Select a single example row and reshape
df_bar <- df %>%
  select(delta_logit_psi, loss_score_orig_fixed) %>%
  slice(1) %>%
  pivot_longer(cols = everything(), names_to = "model", values_to = "score") %>%
  mutate(model = recode(model,
                        delta_logit_psi = "MMSplice (ΔlogitΨ)",
                        loss_score_orig_fixed = "Pangolin (Loss Score)"))

# Step 2: Categorize based on thresholds
df_bar <- df_bar %>%
  mutate(category = case_when(
    model == "MMSplice (ΔlogitΨ)" & score < -2 ~ "high",
    model == "MMSplice (ΔlogitΨ)" & score < -1 ~ "medium",
    model == "MMSplice (ΔlogitΨ)" & score < -0.5 ~ "low",
    model == "Pangolin (Loss Score)" & score < -0.8 ~ "high",
    model == "Pangolin (Loss Score)" & score < -0.5 ~ "medium",
    model == "Pangolin (Loss Score)" & score < -0.2 ~ "low",
    TRUE ~ "none"
  ))

# Step 3: Define colors
category_colors <- c(
  "high" = "red",
  "medium" = "orange",
  "low" = "blue",
  "none" = "grey70"
)

# Step 4: Define cutoff data frames
cutoffs_mmsplice <- data.frame(
  model = "MMSplice (ΔlogitΨ)",
  yintercept = c(-2, -1, -0.5),
  level = c("high", "medium", "low")
)

cutoffs_pangolin <- data.frame(
  model = "Pangolin (Loss Score)",
  yintercept = c(-0.8, -0.5, -0.2),
  level = c("high", "medium", "low")
)

# Step 5: Create plots

# MMSplice
p_mmsplice <- ggplot(df_bar %>% filter(model == "MMSplice (ΔlogitΨ)"),
                     aes(x = model, y = score, fill = category)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_hline(data = cutoffs_mmsplice, aes(yintercept = yintercept, color = level),
             linetype = "dotted", linewidth = 0.8, inherit.aes = FALSE) +
  scale_fill_manual(values = category_colors, name = "Effect") +
  scale_color_manual(values = category_colors, guide = "none") +
  ylim(-2, 0) +
  theme_minimal() +
  theme_cowplot(font_size = fontsize) +
  theme(
    axis.title.x = element_blank(),
    axis.line = element_blank(),
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  labs(
    title = "MMSplice\n(ΔlogitΨ)",
    y = "Score"
    ) +
  theme(
    plot.title = element_text(size = fontsize, face = "plain")
  )
p_mmsplice

# Pangolin
p_pangolin <- ggplot(df_bar %>% filter(model == "Pangolin (Loss Score)"),
                     aes(x = model, y = score, fill = category)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_hline(data = cutoffs_pangolin, aes(yintercept = yintercept, color = level),
             linetype = "dotted", linewidth = 0.8, inherit.aes = FALSE) +
  scale_fill_manual(values = category_colors, name = "Effect") +
  scale_color_manual(values = category_colors, guide = "none") +
  ylim(-1, 0) +
  theme_minimal() +
  theme_cowplot(font_size = fontsize) +
  theme(
    axis.title.x = element_blank(),
    axis.line = element_blank(),
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  labs(
    title = "Pangolin\n(Loss Score)",
    y = "Score"
    ) +
  theme(
    plot.title = element_text(size = fontsize, face = "plain")
  )

p_pangolin  

# ---- Create the legend for the dev example
legend_data <- data.frame(
  variant = "chr1:54723823:C>T",
  junction = "chr1:54723822-54747110:-",
  cutoff = c("high", "medium", "low", "below_low"),
  x = 1:4,
  y = 1
)
legend_data
# First plot: Variant and Junction legends only
legend_plot1 <- ggplot(legend_data, aes(x = x, y = y)) +
  geom_point(aes(shape = variant)) +
  geom_point(aes(fill = junction)) +
  scale_shape_manual(
    name = "Variant",
    values = c("chr1:54723823:C>T" = 16)
  ) +
  scale_linetype_manual(
    name = "Junction",
    values = c("chr1:54723822-54747110:-" = "solid")
  ) +
  theme_void() +
  guides(
    shape = guide_legend(order = 1),
    linetype = guide_legend(order = 2)
  ) +
  theme(legend.position = "right")

legend_plot1

# Second plot: Cutoff legend only
cutoff_colors <- c(high = "red", medium = "orange", low = "blue", below_low = "grey")
# Define a flexible rename map
cutoff_rename <- c(
  high = "high", 
  medium = "medium", 
  low = "low",
  below_low = "below low"
)
legend_data$cutoff_label <- factor(
  cutoff_rename[legend_data$cutoff],
  levels = c("high", "medium", "low", "below low")
)
legend_plot2 <- ggplot(legend_data, aes(x = x, y = y, color = cutoff_label)) +
  geom_point() +
  scale_color_manual(
    name = "Cutoff",
    values = setNames(cutoff_colors, cutoff_rename)
  ) +
  theme_void() +
  guides(color = guide_legend(order = 3)) +
  theme(legend.position = "right")
legend_plot2
# Extract both legends
legend1 <- cowplot::get_legend(legend_plot1)
legend2 <- cowplot::get_legend(legend_plot2)
# Combine the two legends side-by-side
legend_combined <- cowplot::plot_grid(legend1, legend2, ncol = 2, align = "v")
legend_combined


static_row_1row_2col <- ggarrange(
  legend_combined,
  ggarrange(
    p_mmsplice + theme(plot.title = element_text(size = fontsize, face = "plain")),
    p_pangolin + theme(plot.title = element_text(size = fontsize, face = "plain")),
    ncol=2
  ),
  nrow = 1,
  widths = c(2, 1),
  align = "hv"
)
static_row_1row_2col

# ======== STATIC AND DYNAMIC METRICS ===========
final_plot_1col_2col_legend <- ggarrange(
  static_row_1row_2col,
  dev_var_plot,
  ncol = 1,
  heights = c(1, 2.5),
  labels = c('Static Metrics', 'Dynamic Metrics')
)
final_plot_1col_2col_legend


# ================ devAS across LOEUF deciles

df <- fread(file.path(DATA_DIR, 'LOUEF_devAS_distribution.csv'))

devAS_LOEUF_distribution_all_sites <- ggplot(df, aes(x = LOEUF_decile, color = devAS, group = devAS)) +
  geom_point(stat = "count", size = 3, position = position_dodge(width = 0.6)) +
  geom_line(stat = "count", linewidth = 1, position = position_dodge(width = 0.6)) +
  theme_cowplot(font_size = fontsize) +
  labs(
    x = "LOEUF decile", 
    y = "Number of genes",
    color = "Gene contains devAS"
  ) +
  scale_x_continuous(breaks = 1:10) +
  scale_color_manual(values = c("#404040", "#b2182b"))  +
  theme(
    legend.position = c(0.05, 0.05),  # bottom left inside plot
    legend.justification = c(0, 0)
  )
devAS_LOEUF_distribution_all_sites


devAS_LOEUF_distribution_dev_var <- ggplot(df, aes(x = LOEUF_decile, color = dev_var, group = dev_var)) +
  geom_point(stat = "count", size = 3, position = position_dodge(width = 0.6)) +
  geom_line(stat = "count", linewidth = 1, position = position_dodge(width = 0.6)) +
  theme_cowplot(font_size = fontsize) +
  labs(
    x = "LOEUF decile", 
    y = "Number of genes",
    color = "Gene contains Developmental Variant"
  ) +
  scale_x_continuous(breaks = 1:10) +
  scale_color_manual(values = c("#404040", "#b2182b"))  +
  theme(
    legend.position = c(0.2, 0.5),  # bottom left inside plot
    legend.justification = c(0, 0)
  )
devAS_LOEUF_distribution_dev_var



# =============== devAbSplice higher than adult AbSplice

df_num_vars <- fread(file.path(DATA_DIR, 'number_of_variants_high_impact_early_no_impact_adult.csv'))

# Set tissue as an ordered factor by descending num_vars
df_num_vars <- df_num_vars %>%
  filter(tissue != "all_tissues") %>%
  mutate(tissue = factor(tissue, levels = tissue[order(-num_vars)]))

# Create the barplot
devAS_larger_adult <- ggplot(df_num_vars, aes(x = tissue, y = num_vars)) +
  geom_bar(stat = "identity", fill="#b2182b") +
  labs(x = "", y = "Number of Variants") +
  theme_cowplot(font_size = fontsize) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
devAS_larger_adult

df_num_vars <- fread(file.path(DATA_DIR, 'number_of_variants_high_impact_early_no_impact_adult.csv'))
df_num_vars <- df_num_vars %>%
  filter(tissue != "all_tissues") %>%
  mutate(tissue = factor(tissue, levels = tissue[order(num_vars)]))

devAS_larger_adult2 <- ggplot(df_num_vars, aes(y = tissue, x = num_vars)) +
  geom_bar(stat = "identity", fill="#b2182b") +
  labs(y = "", x = "Number of Developmental Variants (high early, low adult)") +
  theme_cowplot(font_size = fontsize)
devAS_larger_adult2


placeholder <- ggplot() + 
  theme_void() +
  theme(panel.background = element_rect(fill = NA, color = NA)) 

# ========== Assemble Fig. 3

fig_3_v4 <- ggarrange(
  ggarrange(
    # dev_example_fig,
    placeholder,
    final_plot_1col_2col_legend,
    nrow=1, labels=c('a', 'b')
  ),
  ggarrange(
    devAS_larger_adult2,
    devAS_LOEUF_distribution_dev_var,
    nrow=1, labels=c('c', 'd')
  ), nrow=2, heights=c(2,0.9)
)
fig_3_v4

# ggsave(
#   filename = file.path(DATA_DIR, "AbSplice2_figures/Fig_3_3.svg"),
#   plot = fig_3_v4,
#   width = 35,
#   height = 20,
#   units = "cm",
#   device = "svg"
# )


fig_3_v5 <- ggarrange(
  ggarrange(
    # dev_example_fig,
    placeholder,
    devAS_larger_adult2,
    devAS_LOEUF_distribution_all_sites,
    ncol=1, labels=c('a', 'c', 'd'), heights=c(4,1.3, 1.3)
    # heights=c(1,1)
  ),
  ggarrange(
    final_plot_1col_2col_legend,
    nrow=1, labels=c('b')
  ), ncol=2 
)
fig_3_v5


ggsave(
  filename = file.path(DATA_DIR, "AbSplice2_figures/Fig_3.svg"),
  plot = fig_3_v5,
  width = 35,
  height = 28,
  units = "cm",
  device = "svg"
)
