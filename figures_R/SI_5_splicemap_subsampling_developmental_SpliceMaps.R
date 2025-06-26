source("config.R")

# ============= SpliceMap subsampling within GTEx =============

results5 <- fread(file.path(DATA_DIR, 'splicemaps_subsampled_liver_gtex.csv'))

library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)

# Define your x values (sample sizes)
x_vals <- c(2, 4, 8, 16, 32, 64)

# Define column renaming
rename_dict <- c(
  ref_psi_correlation = "R^2 of Reference PSI",
  median_n_correlation = "R^2 of median splice site usage",
  jaccard_sites_psi5 = "Jaccard similarity of splice sites"
)

# Define columns to use
columns <- names(rename_dict)

# Assuming results5 is your input data frame

# Aggregate mean and std
grouped_df <- results5 %>%
  filter(num_samples %in% x_vals) %>%
  group_by(num_samples) %>%
  summarise(across(all_of(columns), list(mean = mean, std = sd), .names = "{.col}_{.fn}"), .groups = "drop")

# Reshape for plotting
plot_df <- grouped_df %>%
  pivot_longer(
    cols = -num_samples,
    names_to = c("metric", "stat"),
    names_pattern = "(.*)_(mean|std)"
  ) %>%
  pivot_wider(
    names_from = stat,
    values_from = value
  ) %>%
  mutate(metric_label = rename_dict[metric])


# Plot using ggplot2 with facets
f1 <- ggplot(plot_df, aes(x = num_samples, y = mean)) +
  geom_errorbar(aes(ymin = mean - std, ymax = mean + std), width = 1, color = "grey40") +
  geom_line() +
  geom_point(size = 2) +
  facet_wrap(~ metric_label, scales = "fixed", nrow = 1) +
  scale_x_continuous(breaks = x_vals) +
  ylim(0, 1.1) +
  labs(x = "# of Samples", y = "Value") +
  theme_cowplot(font_size = fontsize)
f1


# ================ Dev with GTEx brain
df_melted <- fread(file.path(DATA_DIR, 'splicemaps_adult_brain_vs_gtex_brains.csv'))

library(ggplot2)
library(dplyr)
library(forcats)

df_melted$tissue <- gsub("_", " ", df_melted$tissue)

# Filter and sort tissues by median for the selected metric
jaccard_median <- df_melted %>%
  filter(Metric == "R² of median splice site usage") %>%
  group_by(tissue) %>%
  summarise(median_val = median(Value), .groups = "drop") %>%
  arrange(median_val)

# Get sorted tissue order
sorted_tissues <- jaccard_median$tissue

# Convert tissue column to factor with correct order
df_melted$tissue <- factor(df_melted$tissue, levels = sorted_tissues)

# Create the plot
f2 <- ggplot(df_melted, aes(x = tissue, y = Value)) +
  geom_boxplot(fill = "lightblue", color = "black", outlier.alpha = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.3, color = "red") +
  facet_wrap(~ Metric, scales = "fixed", nrow = 3) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme_cowplot(font_size = fontsize) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
  ) +
  labs(x = "Tissue", y = "Value")
f2


# ============== Dev with GTEx all tissues
df_melted <- fread(file.path(DATA_DIR, 'splicemaps_adult_all_kaessmann_vs_all_gtex.csv'))

library(ggplot2)
library(dplyr)
library(forcats)

# Filter out kidney tissue
df_melted <- df_melted %>%
  filter(tissue_kaessmann != "kidney")

df_melted$tissue <- gsub("_", " ", df_melted$tissue)

# Compute median Jaccard similarity per tissue for the selected metric
jaccard_median <- df_melted %>%
  filter(Metric == "R² of median splice site usage") %>%
  group_by(tissue_kaessmann) %>%
  summarise(median_val = median(Value), .groups = "drop") %>%
  arrange(median_val)

# Get tissue order from median similarity
tissue_order <- jaccard_median$tissue_kaessmann

# Set the tissue_kaessmann column as a factor with custom order
df_melted <- df_melted %>%
  mutate(tissue_kaessmann = factor(tissue_kaessmann, levels = tissue_order))

# Create the faceted boxplot
f3 <- ggplot(df_melted, aes(x = tissue_kaessmann, y = Value)) +
  geom_boxplot(fill = "lightblue", color = "black", outlier.alpha = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.6, color = "red") +
  facet_wrap(~ Metric, scales = "fixed", nrow = 3) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_cowplot(font_size = fontsize) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
  ) +
  labs(x = "Tissue", y = "Value")

f3

# ================ PR curve GTEx subsampled SpliceMaps
df_performance <- fread(file.path(DATA_DIR, 'pr_curve_gtex_with_subsampled_splicemaps_4_individuals.csv'))

plot_pr_curve <- function(df, color_list, ylim=0.2, breaks_x=0.2, minor_breaks_x=0.2, breaks_y=0.1, minor_breaks_y=0.02, title='All tissues') {
  g <- (
    ggplot(df, aes(x=recall, y=precision, color=model))
    + geom_step(direction = "vh", linewidth=0.9)
    # + scale_color_manual(values=color_list, guide=guide_legend(reverse = TRUE, nrow=4))
    + scale_color_manual(values = color_list, 
                         labels = c("AbSplice-DNA (SpliceAI -> Pangolin)" = "AbSplice-DNA with binary splice site usage", 
                                    "AbSplice-DNA 2" = "AbSplice-DNA with continous splice site usage"),
                         guide = guide_legend(reverse = TRUE, nrow = 4))
    + theme_cowplot(font_size = fontsize)
    #+ background_grid(major = "xy", minor = "xy")
    + scale_x_continuous(limits = c(-0.01, 1.05), breaks = seq(0, 1, breaks_x), minor_breaks = seq(0, 1, minor_breaks_x))
    + scale_y_continuous(limits = c(-0.01, ylim), breaks = seq(0, 1, breaks_y), minor_breaks = seq(0, 1, minor_breaks_y))
    + theme(
      legend.position = c(0.5, 0.8)
      # legend.justification = c("right", "top")
    )
    + labs(x='Recall', y='Precision', color='', title=title)
    + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    # + theme(
    #   legend.position = "top",
    #   #legend.justification = "center"
    # )
    
    # + theme(legend.position = 'none')
  )
  return(g)
}

df_performance <- df_performance %>%
  mutate(model = case_when(
    model == "AbSplice-DNA" ~ "AbSplice-DNA 1\n(SpliceMaps based on all samples)",
    model == "AbSplice_DNA_subsampled_1" ~ "AbSplice-DNA 1\n(SpliceMaps based on 4 samples)",
    TRUE ~ model  # keep other values unchanged
  ))


colors_here <- c(
  "AbSplice-DNA 1\n(SpliceMaps based on all samples)" = "#2ca02c",  # Green
  "AbSplice-DNA 1\n(SpliceMaps based on 4 samples)" = "#d62728"  # Strong Red
)
chosen_models <- c(
  "AbSplice-DNA 1\n(SpliceMaps based on all samples)",
  "AbSplice-DNA 1\n(SpliceMaps based on 4 samples)"
)

df_here <- df_performance[model %in% chosen_models]
df_here <- process_models_pr_curve(df_here, chosen_models)
df_here$model <- factor(df_here$model, levels = chosen_models)

color_list = colors_here[unique(df_here$model)]

pr_curve <- plot_pr_curve(
  df_here,
  ylim=0.40,
  color_list,
  title='Aberrant splicing prediction in GTEx (all tissues)',
)
pr_curve



# ========== All together
fig <- ggarrange(
  ggarrange(f1, pr_curve, ncol=2, labels=c('a', 'b')),
  ggarrange(f2, f3, ncol=2, labels=c('c', 'd')), 
  nrow=2,  heights = c(1,2.5)
)
fig

ggsave(
  filename = file.path(DATA_DIR, "AbSplice2_figures/SI_6_splicemap_subsampling_developmental_SpliceMaps.svg"),
  plot = fig,
  width = 50,
  height = 30,
  units = "cm",
  device = "svg"
)


