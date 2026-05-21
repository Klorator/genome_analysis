# RNA differential expression analysis for the genome assembly project

# Libraries
library(DESeq2)
library(ggplot2)
library(ComplexHeatmap)
library(readr)
library(stringr)
library(dplyr)
library(tibble)
library(rstatix)
library(ggpubr)

# Load data and assemble dataframe
dataFolder <- "2026_VT\\Genome_analysis\\genome_analysis\\result\\07_read_counts"

## From HTSeq-count (canu)
canu_files <- grep("canu", list.files(dataFolder), value = TRUE)

# extract only the condition (bh or serum)
canu_condition <- str_extract(canu_files, "(?<=^rna_).{2,6}(?=_)")
# extract only the sample name part
canu_sample <- str_extract(canu_files, "(?<=^rna_.{2,6}_).+(?=_vs)")

# assemble dataframe with metadata
canu_sampleTable <- data.frame(
  sampleName = canu_sample,
  fileName = canu_files,
  condition = canu_condition
)

canu_sampleTable$condition <- factor(
  canu_sampleTable$condition,
  levels = c("bh", "serum")
)

ddsHTSeq_canu <- DESeqDataSetFromHTSeqCount(
  sampleTable = canu_sampleTable,
  directory = dataFolder,
  design = ~condition
)

# Analysis & results
## canu
ddsHTSeq_canu_results <- DESeq(ddsHTSeq_canu)
res_canu <- results(ddsHTSeq_canu_results)

resultsNames(ddsHTSeq_canu_results)

# log2 FC
resLFC <- lfcShrink(
  ddsHTSeq_canu_results,
  coef = "condition_serum_vs_bh",
  type = "apeglm"
)


# Volcano plot?
plotMA(resLFC, ylim = c(-10, 10)) # kinda crap

logFC_df <- data.frame(
  ID = resLFC@rownames,
  Log2FC = resLFC$log2FoldChange,
  p_val = resLFC$padj
)

cutoff_FC <- 1
cutoff_pval <- 0.05
# add labels for plot fill
logFC_df <- logFC_df |>
  mutate(
    signif = case_when(
      Log2FC > cutoff_FC & p_val < cutoff_pval ~ "serum_upregulated",
      Log2FC < -cutoff_FC & p_val < cutoff_pval ~ "serum_downregulated",
      TRUE ~ "NS"
    ),
    p_val_log = -log10(p_val)
  )

palette <- c(
  "serum_upregulated" = "#E64B35",
  "serum_downregulated" = "#4DBBD5",
  "NS" = "#CCCCCC"
)

## Build volcano plot with ggplot2
vol_plot <- ggplot(logFC_df) +
  aes(x = Log2FC, y = p_val_log) +
  # NS points
  geom_point(
    data = subset(logFC_df, signif == "NS"),
    color = palette["NS"],
    alpha = 0.3,
    show.legend = FALSE
  ) +
  # Reference lines for significance thresholds
  geom_vline(
    xintercept = c(-cutoff_FC, cutoff_FC),
    color = palette["NS"],
    alpha = 0.5
  ) +
  geom_hline(
    yintercept = -log10(cutoff_pval),
    color = palette["NS"],
    alpha = 0.5
  ) +
  # Upregulated points
  geom_point(
    data = subset(logFC_df, signif == "serum_upregulated"),
    aes(color = signif),
    # color = palette["serum_upregulated"],
    alpha = 0.7
  ) +
  # Downregulated points
  geom_point(
    data = subset(logFC_df, signif == "serum_downregulated"),
    aes(color = signif),
    # color = palette["serum_downregulated"],
    alpha = 0.7
  ) +
  scale_color_manual(
    values = c(
      serum_upregulated = "#E64B35",
      serum_downregulated = "#4DBBD5"
    ),
    breaks = c("serum_upregulated", "serum_downregulated"),
    labels = c("Serum", "BH"),
    name = NULL,
    drop = FALSE
  ) +
  # Axis limits
  xlim(c(-10, 10)) +
  ylim(c(0, 325)) + # values above 325 are INF due to p_val = 0
  # Labels and theme
  labs(
    title = NULL,
    x = "Log2 Fold Change (BH vs Serum)",
    y = "-Log10 Adjusted P-value"
  ) +
  theme_classic() +
  # Increase axis labels and ticks size
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

ggsave(
  filename = "2026_VT\\Genome_analysis\\genome_analysis\\result\\volcano_plot.png",
  plot = vol_plot,
  width = 8,
  height = 6,
  dpi = 300
)

###################################
# Heatmap
heatmap_df <- counts(ddsHTSeq_canu_results, normalized = TRUE) |>
  as.data.frame() |>
  rownames_to_column("gene_id") |>
  na.omit()
# Zscore normalization for better visualization
heatmap_df[, -1] <- t(apply(heatmap_df[, -1], 1, scale))
heatmap_df <- heatmap_df[
  apply(heatmap_df[, -1], 1, function(x) all(is.finite(x))),
  ,
  drop = FALSE
]

sample_anno <- HeatmapAnnotation(
  condition = canu_sampleTable$condition,
  col = list(condition = c("serum" = "#E64B35", "bh" = "#4DBBD5")),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    condition = list(
      title = "Condition",
      at = c("serum", "bh"),
      labels = c("Serum", "BH")
    )
  )
)


heatmap_full <- Heatmap(
  as.matrix(heatmap_df[, -1]),
  name = "Normalized Counts",
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_names_rot = 45,
  width = unit(6, "in"),
  top_annotation = sample_anno,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_dend = FALSE,
  col = circlize::colorRamp2(c(-2, 0, 2), c("#00204C", "#7E7F7A", "#FDE333")),
  heatmap_legend_param = list(
    title = "Normalized counts\n(Z-score)",
    at = c(-2, 0, 2)
  )
)

png(
  filename = "2026_VT\\Genome_analysis\\genome_analysis\\result\\heatmap.png",
  width = 8,
  height = 10,
  units = "in",
  res = 300
)
draw(
  heatmap_full,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legend = TRUE
)
dev.off()

##################################################
# Barplot of genes of interest with stats
## Genes of interest
genes_of_interest <- data.frame(
  Gene = c(
    "manZ_3",
    "manY_2",
    "ptsL",
    "algB",
    "guaB",
    "purA",
    "pyrF",
    "pyrK_2",
    "purD",
    "purH",
    "purL",
    "purQ",
    "purC",
    "clsA_1",
    "ddcP",
    "ldt_fmp",
    "mgs",
    "lytA_2"
  ),
  Description = c(
    "Mannose PTS permease subunit (sugar transport)",
    "Mannose PTS permease subunit (sugar transport)",
    "Phosphotransferase system (PTS) component (sugar uptake)",
    "Response regulator / transcriptional regulator",
    "IMP dehydrogenase (purine biosynthesis)",
    "Adenylosuccinate synthase (purine biosynthesis)",
    "Orotidine 5'-phosphate decarboxylase (pyrimidine biosynthesis)",
    "Dihydroorotate dehydrogenase-related subunit (pyrimidine biosynthesis)",
    "Phosphoribosylamine–glycine ligase (purine biosynthesis)",
    "AICAR transformylase/IMP cyclohydrolase (purine biosynthesis)",
    "Phosphoribosylformylglycinamidine synthase (purine biosynthesis)",
    "Phosphoribosylformylglycinamidine synthase subunit (purine biosynthesis)",
    "SAICAR synthase (purine biosynthesis)",
    "Cardiolipin synthase (membrane lipid biosynthesis)",
    "D-alanyl-D-alanine carboxypeptidase / cell wall remodeling",
    "L,D-transpeptidase (peptidoglycan cross-linking)",
    "Methylglyoxal synthase (glycolytic side-path)",
    "Autolysin (peptidoglycan hydrolase)"
  ),
  Group = c(
    "PTS / Sugar Transport",
    "PTS / Sugar Transport",
    "PTS / Sugar Transport",
    "Regulation",
    "Purine Biosynthesis I",
    "Purine Biosynthesis I",
    "Pyrimidine Biosynthesis",
    "Pyrimidine Biosynthesis",
    "Purine Biosynthesis I",
    "Purine Biosynthesis I",
    "Purine Biosynthesis II",
    "Purine Biosynthesis II",
    "Purine Biosynthesis II",
    "Cell Wall / Envelope",
    "Cell Wall / Envelope",
    "Cell Wall / Envelope",
    "PTS / Sugar Transport",
    "Cell Wall / Envelope"
  ),
  stringsAsFactors = FALSE
)
gff <- rtracklayer::import(file.path(
  dataFolder,
  grep(".gff", list.files(dataFolder), value = TRUE)
))
gff_df <- as.data.frame(mcols(gff)) |> as_tibble()

# Join gff annotation with logFC_df to get gene names
logFC_df_annotated <- logFC_df |>
  left_join(gff_df[, c("locus_tag", "gene")], by = c("ID" = "locus_tag"))
# Join the normalized counts from heatmap_df
sample_cols <- colnames(heatmap_df)[-1]
logFC_df_annotated <- logFC_df_annotated |>
  left_join(
    heatmap_df[, c("gene_id", sample_cols)],
    by = c("ID" = "gene_id")
  )
# Filter for genes of interest
logFC_genes_of_interest <- logFC_df_annotated |>
  filter(gene %in% genes_of_interest$Gene) |>
  left_join(genes_of_interest, by = c("gene" = "Gene"))
# Cleanup duplicated columns and reorder
logFC_genes_of_interest <- logFC_genes_of_interest |>
  mutate(
    Log2FC = Log2FC.x,
    p_val = p_val.x,
    signif = signif.x
  ) |>
  select(
    gene,
    Description,
    Group,
    Log2FC,
    p_val,
    signif,
    any_of(sample_cols)
  ) |>
  arrange(Group, gene)

# Average and standard deviation columns for each sample set
serum_samples <- canu_sampleTable |>
  filter(condition == "serum") |>
  pull(sampleName)

bh_samples <- canu_sampleTable |>
  filter(condition == "bh") |>
  pull(sampleName)

logFC_genes_of_interest <- logFC_genes_of_interest |>
  mutate(
    serum_avg = rowMeans(pick(all_of(serum_samples)), na.rm = TRUE),
    bh_avg = rowMeans(pick(all_of(bh_samples)), na.rm = TRUE),
    serum_sd = apply(pick(all_of(serum_samples)), 1, sd, na.rm = TRUE),
    bh_sd = apply(pick(all_of(bh_samples)), 1, sd, na.rm = TRUE)
  )

# Stars for significance
logFC_genes_of_interest <- logFC_genes_of_interest |>
  mutate(
    stars = case_when(
      signif == "NS" ~ NA,
      p_val < 0.001 ~ "***",
      p_val < 0.01 ~ "**",
      p_val < 0.05 ~ "*",
      TRUE ~ NA
    )
  )

# Create dataframe for stats brackets using rstatix
stats_df <- logFC_genes_of_interest |>
  select(gene, stars, any_of(sample_cols)) |>
  tidyr::pivot_longer(
    cols = any_of(sample_cols),
    names_to = "sampleName",
    values_to = "avg"
  ) |>
  left_join(
    canu_sampleTable |> select(sampleName, condition),
    by = "sampleName"
  ) |>
  mutate(condition = factor(condition, levels = c("serum", "bh"))) |>
  group_by(gene) |>
  t_test(avg ~ condition) |>
  add_xy_position(x = "gene", group = "condition", dodge = 0.8) |>
  left_join(
    logFC_genes_of_interest |> select(gene, stars),
    by = "gene"
  ) |>
  mutate(
    label = stars,
    p = NA_real_,
    p.adj = NA_real_,
    p.signif = stars,
    p.adj.signif = stars
  )

# Prepare data for barplot using normalized counts
anno <- logFC_df_annotated |>
  select(
    ID,
    gene,
    dplyr::any_of(c("Group", "group", "GroupName", "group_name"))
  ) |>
  distinct()

if (!"Group" %in% colnames(anno)) {
  if ("group" %in% colnames(anno)) {
    anno <- anno |> rename(Group = group)
  } else if ("GroupName" %in% colnames(anno)) {
    anno <- anno |> rename(Group = GroupName)
  } else if ("group_name" %in% colnames(anno)) {
    anno <- anno |> rename(Group = group_name)
  } else {
    anno$Group <- NA_character_
  }
}

normalized_counts_long <- counts(ddsHTSeq_canu_results, normalized = TRUE) |>
  as.data.frame() |>
  rownames_to_column("ID") |>
  tidyr::pivot_longer(
    cols = -ID,
    names_to = "sampleName",
    values_to = "count"
  ) |>
  left_join(anno, by = "ID") |>
  filter(gene %in% logFC_genes_of_interest$gene) |>
  left_join(
    canu_sampleTable |> select(sampleName, condition),
    by = "sampleName"
  ) |>
  mutate(condition = factor(condition, levels = c("serum", "bh")))

# Calculate summary statistics per gene and condition
barplot_summary <- normalized_counts_long |>
  group_by(gene, Group, condition) |>
  summarise(
    avg = mean(count, na.rm = TRUE),
    sd = sd(count, na.rm = TRUE),
    .groups = "drop"
  )

# Combine with individual samples for plotting
barplot_data <- normalized_counts_long |>
  rename(sample_count = count) |>
  left_join(
    barplot_summary,
    by = c("gene", "Group", "condition")
  ) |>
  select(gene, Group, condition, sampleName, sample_count, avg, sd)

# Build barplot
barplot <- ggplot(barplot_data, aes(x = gene, y = avg, fill = condition)) +
  geom_col(position = position_dodge(0.8), alpha = 0.8) +
  geom_errorbar(
    aes(ymin = avg - sd, ymax = avg + sd),
    position = position_dodge(0.8),
    width = 0.2,
    linewidth = 0.5
  ) +
  geom_point(
    aes(y = sample_count, color = condition),
    position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.22),
    alpha = 0.5,
    size = 3,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = c("serum" = "#E64B35", "bh" = "#4DBBD5"),
    labels = c("Serum", "BH")
  ) +
  scale_color_manual(
    values = c("serum" = "#B53A28", "bh" = "#2F8CA8")
  ) +
  scale_y_log10(labels = scales::label_number()) +
  annotation_logticks(sides = "l") +
  labs(
    title = NULL,
    x = "Gene",
    y = "Normalized Counts (log10)",
    fill = "Condition"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )

barplot

ggsave(
  filename = "2026_VT\\Genome_analysis\\genome_analysis\\result\\barplot_genes_of_interest.png",
  plot = barplot,
  width = 14,
  height = 8,
  dpi = 300
)
