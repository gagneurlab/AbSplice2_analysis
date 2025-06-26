source("config.R")

# ===================== Across gene TPM =====================

colors_here <- c(
  "SpliceAI" = "#6b9dd4",         # Medium Blue (clearer than #b3d9ff)
  "Pangolin" = "#e69f00",         # Warm Orange (improves contrast from yellow)
  "AbSplice-DNA 1" = "#e41a1c",   # Muted Pink/Magenta (distinct but not too bright)
  "AbSplice-DNA 2" = "#a50f15",   # Strong Red (keeping it distinct)
  "AbSplice-DNA with gene TPM" = "#984ea3"  # Medium Purple (softer, less harsh)
)

chosen_models <- c(
  "SpliceAI",
  "Pangolin",
  "AbSplice-DNA 1",
  "AbSplice-DNA 2",
  "AbSplice-DNA with gene TPM"
)

df <- fread(pr_curve_across_gene_expression_values)
df <- process_models_pr_curve(df, chosen_models)
df_stats <- fread(pr_curve_across_gene_expression_values_stats)
color_list = colors_here[unique(df$model)]

# define the order of the facets
tpm_cutoff_order <- c(
  "(0, 3)",    
  "(3, 5)",    
  "(5, 10)",   
  "(10, 20)",  
  "(20, 50)",  
  "(50, inf)"
)
df_stats$tpm_cutoff <- factor(df_stats$tpm_cutoff, levels=tpm_cutoff_order)
df$cutoff_pair <- factor(df$cutoff_pair, levels = tpm_cutoff_order)
df_mapping <- df_stats %>% select(tpm_cutoff, tpm_cutoff_label) %>% unique()
mapping_list <- setNames(df_mapping$tpm_cutoff_label, df_mapping$tpm_cutoff)

label_dict <- c(
  "SpliceAI" = "SpliceAI",
  "Pangolin" = "Pangolin",
  "AbSplice-DNA 1" = "AbSplice1",
  "AbSplice-DNA 2" = "AbSplice2",
  "AbSplice-DNA with gene TPM" = "AbSplice with gene TPM"
)

df$model_label <- factor(label_dict[as.character(df$model)], levels = label_dict[chosen_models])
color_list <- setNames(colors_here[names(label_dict)], label_dict[names(label_dict)])


plot_pr_curve_variant <- function(df, color_list, ylim=0.2, breaks_x=0.2, minor_breaks_x=0.2, breaks_y=0.1, minor_breaks_y=0.02, title='All tissues') {
  g <- (
    ggplot(df, aes(x=recall, y=precision, color=model_label))
    + geom_step(direction = "vh", linewidth=0.9)
    + scale_color_manual(
      values=color_list, 
      guide=guide_legend(
        reverse = TRUE, nrow=2
      )
    )
    # + scale_color_manual(values=color_list, guide=guide_legend(reverse = TRUE, nrow=4))
    + theme_cowplot(font_size = fontsize)
    # + background_grid(major = "xy", minor = "xy")
    + scale_x_continuous(limits = c(-0.01, 1.05), breaks = seq(0, 1, breaks_x), minor_breaks = seq(0, 1, minor_breaks_x))
    + scale_y_continuous(limits = c(-0.01, ylim), breaks = seq(0, 1, breaks_y), minor_breaks = seq(0, 1, minor_breaks_y))
    + labs(x='Recall', y='Precision', color='', title=title)
    + facet_wrap('cutoff_pair', nrow=2, scales='free_y', labeller=labeller(cutoff_pair=mapping_list))
    # + theme(legend.text = element_text(
    #   margin = margin(r = 30, unit = "pt")))
    # + panel_border()
    + theme(
      legend.key = element_blank(),
      legend.position = "bottom",
      strip.background = element_rect(colour="white", fill="white")
    )
  )
  return(g)
}

fig <- plot_pr_curve_variant(
  df,
  ylim=1.0,
  color_list,
  title=''
)
fig

ggsave(
  filename = file.path(DATA_DIR, "AbSplice2_figures/SI_4_across_tpm_values.svg"),
  plot = fig,
  width = 30,
  height = 25,
  units = "cm",
  device = "svg"
)
