library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

source("config.R")

# Load data
df <- fread(file.path(DATA_DIR, 'developmental_variant_example__discontinuity.csv'))
chosen_variant <- "chr20:54948463:C>G"
chosen_junction <- "chr20:54945715-54948463:-"


df <- fread(file.path(DATA_DIR, 'developmental_variant_example__discontinuity2.csv'))
chosen_variant <- 'chr4:493281:A>C'
chosen_junction <- 'chr4:493278-494184:+'

# Define valid timepoints and their labels
timepoints <- c(paste0("t", 1:15), "tgtex")
timepoint_labels <- c(
  t1 = "4 wpc", t2 = "5 wpc", t3 = "6 wpc", t4 = "7 wpc", t5 = "8 wpc",
  t6 = "9–11 wpc", t7 = "12 wpc", t8 = "13 wpc", t9 = "16 wpc",
  t10 = "18–19 wpc", t11 = "20 wpc", t12 = "Newborn–Toddler",
  t13 = "School–Young Adult", t14 = "Young Mid Age", t15 = "Older Mid–Senior",
  tgtex = "GTEx"
)

# Filter out unexpected values
df <- df %>% filter(timepoint %in% timepoints)

# Assign factor and relabel timepoints
df$timepoint <- factor(df$timepoint, levels = timepoints, labels = timepoint_labels)

# Create cutoff-based categories for both AbSplice1 and AbSplice2
df <- df %>%
  mutate(
    absplice2_category = case_when(
      AbSplice2 > 0.2 ~ "high",
      AbSplice2 > 0.1 ~ "medium",
      AbSplice2 > 0.05 ~ "low",
      TRUE ~ "none"
    ),
    absplice1_category = case_when(
      AbSplice1 > 0.2 ~ "high",
      AbSplice1 > 0.1 ~ "medium",
      AbSplice1 > 0.05 ~ "low",
      TRUE ~ "none"
    ),
    delta_psi_category = case_when(
      delta_psi <= -0.2 ~ "high",
      delta_psi <= -0.1 ~ "medium",
      delta_psi <= -0.05 ~ "low",
      TRUE ~ "none"
    )
  )

# Reshape to long format with all metrics of interest
df_long <- df %>%
  pivot_longer(
    cols = c(ref_psi, AbSplice1, AbSplice2, delta_psi, median_n, splice_site_is_expressed),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    absplice_category = case_when(
      metric == "AbSplice2" ~ absplice2_category,
      metric == "AbSplice1" ~ absplice1_category,
      metric == "delta_psi" ~ delta_psi_category,
      TRUE ~ NA_character_
    ),
    facet_metric = recode(metric,
                          ref_psi = "Reference PSI",
                          AbSplice1 = "AbSplice1",
                          AbSplice2 = "AbSplice2",
                          delta_psi = "ΔPSI",
                          median_n = "Median # Splice Site Reads",
                          splice_site_is_expressed = "Splice Site Expressed"
    ),
    facet_metric = factor(
      facet_metric,
      levels = c(
        "Reference PSI", "AbSplice1", "AbSplice2",
        "ΔPSI", "Median # Splice Site Reads", "Splice Site Expressed"
      )
    )
  )

# Define color mapping
category_colors <- c(
  "high" = "red",
  "medium" = "orange",
  "low" = "blue",
  "none" = "grey70"
)


# Horizontal cutoff lines for AbSplice1, AbSplice2, and ΔPSI
cutoffs_df <- data.frame(
  facet_metric = rep(c("AbSplice1", "AbSplice2", "ΔPSI"), each = 3),
  yintercept = c(0.05, 0.1, 0.2, 0.05, 0.1, 0.2, -0.05, -0.1, -0.2),
  cutoff_level = factor(rep(c("low", "medium", "high"), 3), levels = c("high", "medium", "low")),
  color = rep(c("blue", "orange", "red"), 3)
)
cutoffs_df$facet_metric <- factor(cutoffs_df$facet_metric, levels = levels(df_long$facet_metric))





# Fix y-axis range for Reference PSI
psi_limits_df <- data.frame(
  timepoint = factor("4 wpc", levels = timepoint_labels),
  value = c(0, 1),
  facet_metric = factor("Reference PSI", levels = levels(df_long$facet_metric))
)
# Force ΔPSI y-axis range
delta_psi_range <- range(df$delta_psi, na.rm = TRUE)
delta_psi_limits_df <- data.frame(
  timepoint = factor("4 wpc", levels = timepoint_labels),
  value = c(min(-0.2, delta_psi_range[1]), max(0.2, delta_psi_range[2])),
  facet_metric = factor("ΔPSI", levels = levels(df_long$facet_metric))
)


metric_order <- c(
  "Reference PSI",
  "Splice Site Expressed",
  "Median # Splice Site Reads",
  "ΔPSI",
  "AbSplice1",
  "AbSplice2"
)

df_long <- df_long %>%
  mutate(
    facet_metric = factor(facet_metric, levels = metric_order)
  )


cutoffs_mediann_df <- data.frame(
  timepoint = factor("4 wpc", levels = timepoint_labels),
  yintercept = 10,
  facet_metric = factor("Median # Splice Site Reads", levels = levels(df_long$facet_metric))
)


# Final plot
dev_var_plot <- ggplot(df_long, aes(x = timepoint, y = value)) +
  # Force y-axis range for Reference PSI
  geom_blank(data = psi_limits_df, aes(x = timepoint, y = value)) +
  geom_blank(data = delta_psi_limits_df, aes(x = timepoint, y = value)) +
  
  # Colored points for AbSplice1 and AbSplice2
  geom_point(
    data = df_long %>% filter(metric %in% c("AbSplice1", "AbSplice2", "delta_psi")),
    aes(color = absplice_category),
    size = 2,
    show.legend = FALSE
  ) +
  
  # Horizontal cutoff lines
  geom_hline(
    data = cutoffs_df,
    aes(yintercept = yintercept, color = cutoff_level),
    linetype = "dotted",
    show.legend = FALSE
  ) +
  
  geom_hline(
    data = cutoffs_mediann_df,
    aes(yintercept = yintercept),
    color = "grey50",
    linetype = "dashed",
    linewidth = 0.6,
    inherit.aes = FALSE
  ) +
  
  
  # Black points for all other metrics
  geom_point(
    data = df_long %>% filter(!(metric %in% c("AbSplice1", "AbSplice2", "delta_psi"))),
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
    title = "",
    y = "Value",
    x = "Timepoint"
  ) +
  
  facet_wrap(vars(facet_metric), scales = "free_y", nrow = 2)

# Show the plot
dev_var_plot



# ============ STATIC ============

# ================================================  Static scores
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

# Step 1: Select a single example row and reshape
df_bar <- df %>%
  select(delta_logit_psi, loss_score, delta_score) %>%
  slice(1) %>%
  pivot_longer(cols = everything(), names_to = "model", values_to = "score") %>%
  mutate(model = recode(model,
                        delta_logit_psi = "MMSplice (ΔlogitΨ)",
                        loss_score = "Pangolin (Loss Score)",
                        delta_score = "SpliceAI (Delta Score)"))

df_bar
# Step 2: Categorize based on thresholds
df_bar <- df_bar %>%
  mutate(category = case_when(
    model == "MMSplice (ΔlogitΨ)" & score < -2 ~ "high",
    model == "MMSplice (ΔlogitΨ)" & score < -1 ~ "medium",
    model == "MMSplice (ΔlogitΨ)" & score < -0.5 ~ "low",
    model == "Pangolin (Loss Score)" & score < -0.8 ~ "high",
    model == "Pangolin (Loss Score)" & score < -0.5 ~ "medium",
    model == "Pangolin (Loss Score)" & score < -0.2 ~ "low",
    model == "SpliceAI (Delta Score)" & score > 0.8 ~ "high",
    model == "SpliceAI (Delta Score)" & score > 0.5 ~ "medium",
    model == "SpliceAI (Delta Score)" & score > 0.2 ~ "low",
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

cutoffs_spliceai <- data.frame(
  model = "SpliceAI (Delta Score)",
  yintercept = c(0.8, 0.5, 0.2),
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
  # coord_flip() +
  ylim(-5, 0) +
  theme_minimal() +
  theme_cowplot(font_size = fontsize) +
  theme(
    axis.title.x = element_blank(),
    # axis.text.y = element_text(face = "bold"),
    axis.line = element_blank(),
    # legend.position = "right"
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  # theme(
  #   axis.title.y = element_blank(),
  #   axis.text.y = element_text(face = "bold"),
  #   legend.position = "none"
  # ) +
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
  # coord_flip() +
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


# SpliceAI
p_spliceai <- ggplot(df_bar %>% filter(model == "SpliceAI (Delta Score)"),
                     aes(x = model, y = score, fill = category)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_hline(data = cutoffs_spliceai, aes(yintercept = yintercept, color = level),
             linetype = "dotted", linewidth = 0.8, inherit.aes = FALSE) +
  scale_fill_manual(values = category_colors, name = "Effect") +
  scale_color_manual(values = category_colors, guide = "none") +
  # coord_flip() +
  ylim(0, 1) +
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
    title = "SpliceAI\n(Delta Score)",
    y = "Score"
  ) +
  theme(
    plot.title = element_text(size = fontsize, face = "plain")
  )

p_spliceai 


library(ggplot2)
library(cowplot)

# Create dummy data if needed
legend_data <- data.frame(
  variant = chosen_variant,
  junction = chosen_junction,
  cutoff = c("high", "medium", "low", "below_low"),
  x = 1:4,
  y = 1
)
legend_data

legend_plot1 <- ggplot(legend_data, aes(x = x, y = y)) +
  geom_point(aes(shape = variant), size = 2) +
  geom_point(aes(fill = junction), shape = 21, size = 2, color = "black") +
  scale_shape_manual(
    name = "Variant",
    values = setNames(16, chosen_variant)
  ) +
  scale_fill_manual(
    name = "Junction",
    values = setNames("black", chosen_junction)
  ) +
  theme_void() +
  guides(
    shape = guide_legend(order = 1),
    fill = guide_legend(order = 2)
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

# Apply rename directly to the factor levels
# legend_data$cutoff_label <- cutoff_rename[legend_data$cutoff]
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
  # ggdraw() + draw_grob(legend_combined, x = 0, y = 0, width = 1, height = 1),
  legend_combined,
  ggarrange(
    p_mmsplice + theme(plot.title = element_text(size = fontsize, face = "plain")),
    p_spliceai + theme(plot.title = element_text(size = fontsize, face = "plain")),
    p_pangolin + theme(plot.title = element_text(size = fontsize, face = "plain")),
    ncol=3
  ),
  nrow = 1,
  widths = c(1, 1),
  align = "hv"
)
static_row_1row_2col


# ======== STATIC AND DYNAMIC METRICS ===========
final_plot_1col_2col_legend <- ggarrange(
  # static_with_title,
  static_row_1row_2col,
  dev_var_plot,
  ncol = 1,
  heights = c(1, 2.5),
  labels = c('Static Metrics', 'Dynamic Metrics')
)
final_plot_1col_2col_legend


ggsave(
  filename = file.path(DATA_DIR, "AbSplice2_figures/SI_9_developmental_variant__discontinuity.svg"),
  plot = final_plot_1col_2col_legend,
  width = 30,
  height = 20,
  units = "cm",
  device = "svg"
)




