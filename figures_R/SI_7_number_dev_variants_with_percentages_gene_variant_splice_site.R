library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)

# Load data
df_combined <- fread(file.path(DATA_DIR, 'developmental_predictions_stats_plot.csv'))

# Prepare
df_combined <- df_combined %>%
  mutate(
    type = factor(type, levels = c("Purely Developmental", "Developmental", "Adult", "Purely Adult", "Total")),
    class = factor(class, levels = c("Variants", "Genes", "Splice Sites"))
  )

# Sort tissues based on 'Developmental' variant counts
sort_order <- df_combined %>%
  filter(type == "Developmental", class == "Variants") %>%
  arrange(desc(count)) %>%
  pull(group) %>% unique()

df_combined$group <- factor(df_combined$group, levels = sort_order)

# Plot
fig_number_variants <- ggplot(df_combined, aes(x = group, y = count)) +
  geom_bar(
    data = df_combined %>% filter(type != "Total"),
    aes(fill = type),
    stat = "identity",
    position = position_dodge(width = 0.9)
  ) +
  geom_line(
    data = df_combined %>% filter(type == "Total"),
    aes(group = type, linetype = type, color = type),
    size = 0.6
  ) +
  geom_point(
    data = df_combined %>% filter(type == "Total"),
    aes(shape = type, color = type),
    size = 2
  ) +
  facet_wrap(~class, scales = "free_x", ncol = 1) +
  coord_flip() +
  scale_fill_manual(
    name = "Type",
    values = c(
      "Purely Developmental" = "#5a3d84",
      "Developmental" = "#8c6bb1",
      "Adult" = "#fdb863",
      "Purely Adult" = "#e08214"
    )
  ) +
  scale_color_manual(
    name = "Type",
    values = c("Total" = "black"),
    breaks = "Total"
  ) +
  scale_linetype_manual(
    name = "Type",
    values = c("Total" = "dotted"),
    breaks = "Total"
  ) +
  scale_shape_manual(
    name = "Type",
    values = c("Total" = 16),
    breaks = "Total"
  ) +
  labs(x = NULL, y = "Count", title = "") +
  theme_cowplot() +
  theme(
    axis.text.y = element_text(size = 10),
    legend.position = "right"
  )

fig_number_variants

ggsave(
  filename = file.path(DATA_DIR, "AbSplice2_figures/SI_Fig_number_developmental_variants_genes_variants_splice_sites.svg"),
  plot = fig_number_variants,
  width = 25,
  height = 25,
  units = "cm",
  device = "svg"
)



