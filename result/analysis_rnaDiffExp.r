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
    x = "Log2 Fold Change (Serum vs BH)",
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
  col = circlize::colorRamp2(c(-2, 0, 2), hcl_palette = "Cividis"),
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
genes_of_interest <- data.frame(
  Gene = c(
    # PTS / sugar transport
    "manZ_3",
    "manY_2",
    "ptsL",
    "ptsI",
    # Regulation
    "algB",
    "afr_2",
    "rpoN",
    "IsrC_1",
    # Purine biosynthesis
    "guaB",
    "purA",
    "purD",
    "purH",
    "purL",
    "purQ",
    "purC",
    # Pyrimidine biosynthesis
    "pyrF",
    "pyrK_2",
    # Cell wall / envelope & membrane
    "clsA_1",
    "ddcP",
    "ldt_fmp",
    "mgs",
    "lytA_2",
    # Transport / uptake (non‑PTS)
    "artM_1",
    "bioY2",
    "hmpT",
    # Other Tn‑seq hits (hypothetical / uncharacterized)
    "hypBA2",
    "EfmE745_03139",
    "EfmE745_03220",
    "EfmE745_03141",
    "EfmE745_03101",
    "EfmE745_03147",
    "EfmE745_03131",
    "EfmE745_01958",
    "EfmE745_01972",
    "EfmE745_02302",
    "EfmE745_00881",
    "rpsN2"
  ),

  Description = c(
    # PTS / sugar transport
    "Mannose PTS permease subunit (IIBC component, sugar transport)",
    "Mannose PTS permease subunit (IIA component, sugar transport)",
    "PTS enzyme I component (sugar uptake phosphotransfer)",
    "PTS enzyme I component (canonical ptsI annotation)",
    # Regulation
    "Putative transcriptional regulator (AlgB-like)",
    "Putative transcriptional regulator (AraC-family)",
    "Sigma-54 (RpoN) alternative sigma factor",
    "Putative transcriptional regulator (LysR/other)",
    # Purine biosynthesis
    "IMP dehydrogenase (GMP branch, purine biosynthesis)",
    "Adenylosuccinate synthase (AMP branch, purine biosynthesis)",
    "Phosphoribosylamine–glycine ligase (PurD, de novo purine)",
    "AICAR transformylase/IMP cyclohydrolase (PurH, de novo purine)",
    "FGAM synthase (PurL, de novo purine)",
    "FGAM synthase subunit (PurQ, de novo purine)",
    "SAICAR synthase (PurC, de novo purine)",
    # Pyrimidine biosynthesis
    "Orotidine 5'-phosphate decarboxylase (PyrF, de novo pyrimidine)",
    "Dihydroorotate dehydrogenase-related subunit (PyrK_2, de novo pyrimidine)",
    # Cell wall / envelope & membrane
    "Cardiolipin synthase (membrane phospholipid remodeling)",
    "D-alanyl-D-alanine carboxypeptidase (low-MW PBP, cell wall remodeling)",
    "L,D-transpeptidase (3→3 peptidoglycan cross-linking)",
    "Monoglucosyldiacylglycerol synthase (membrane glycolipid synthesis)",
    "Autolysin (peptidoglycan hydrolase)",
    # Transport / uptake (non‑PTS)
    "ABC transporter permease (arginine/AA or peptide uptake)",
    "Biotin transporter (BioY family)",
    "Transporter (HMP family, small molecule uptake)",
    # Other Tn‑seq hits
    "Carbohydrate metabolism protein (HypBA2-like)",
    "Hypothetical protein (Tn‑seq fitness hit)",
    "Hypothetical protein (Tn‑seq fitness hit)",
    "Hypothetical protein (Tn‑seq fitness hit)",
    "Hypothetical protein (Tn‑seq fitness hit)",
    "Hypothetical protein (Tn‑seq fitness hit)",
    "Hypothetical protein (Tn‑seq fitness hit)",
    "Hypothetical protein (Tn‑seq fitness hit)",
    "Hypothetical protein (Tn‑seq fitness hit)",
    "Hypothetical protein (Tn‑seq fitness hit)",
    "Hypothetical protein (Tn‑seq fitness hit)",
    "Ribosomal protein S14 paralog (rpsN2)"
  ),

  Group = c(
    # PTS / sugar transport
    "PTS / Sugar Transport",
    "PTS / Sugar Transport",
    "PTS / Sugar Transport",
    "PTS / Sugar Transport",
    # Regulation
    "Regulation",
    "Regulation",
    "Regulation",
    "Regulation",
    # Purine biosynthesis
    "Purine Biosynthesis",
    "Purine Biosynthesis",
    "Purine Biosynthesis",
    "Purine Biosynthesis",
    "Purine Biosynthesis",
    "Purine Biosynthesis",
    "Purine Biosynthesis",
    # Pyrimidine biosynthesis
    "Pyrimidine Biosynthesis",
    "Pyrimidine Biosynthesis",
    # Cell wall / envelope & membrane
    "Cell Wall / Envelope",
    "Cell Wall / Envelope",
    "Cell Wall / Envelope",
    "Cell Wall / Envelope",
    "Cell Wall / Envelope",
    # Transport / uptake (non‑PTS)
    "Transport",
    "Transport",
    "Transport",
    # Other Tn‑seq hits
    "Carbohydrate Metabolism",
    "Other",
    "Other",
    "Other",
    "Other",
    "Other",
    "Other",
    "Other",
    "Other",
    "Other",
    "Other",
    "Ribosome"
  ),

  stringsAsFactors = FALSE
)

###################
genes_subsets <- list(
  ## Core metabolic bottlenecks (nucleotide + PTS)
  goi_nucleotide = c(
    "guaB",
    "purA",
    "purD",
    "purH",
    "purL",
    "purQ",
    "purC",
    "pyrF",
    "pyrK_2"
  ),

  goi_pts = c(
    "manZ_3",
    "manY_2",
    "ptsL",
    "ptsI"
  ),

  ## Cell wall / envelope genes (including negative-fitness hits)
  goi_cellwall = c(
    "clsA_1",
    "ddcP",
    "ldt_fmp",
    "mgs",
    "lytA_2"
  ),

  ## Regulators (good for “control layer” plots)
  goi_regulators = c(
    "algB",
    "afr_2",
    "rpoN",
    "IsrC_1"
  ),

  ## Transporters beyond PTS
  goi_transport_nonPTS = c(
    "artM_1",
    "bioY2",
    "hmpT"
  ),

  ## “Negative fitness” genes (mutants grow better in serum)
  goi_negative_fitness = c(
    "clsA_1",
    "ddcP",
    "ldt_fmp",
    "mgs",
    "lytA_2"
  ),

  ## “Positive fitness” genes (required for growth in serum)
  goi_positive_fitness = c(
    "manZ_3",
    "manY_2",
    "ptsL",
    "ptsI",
    "algB",
    "afr_2",
    "rpoN",
    "IsrC_1",
    "guaB",
    "purA",
    "purD",
    "purH",
    "purL",
    "purQ",
    "purC",
    "pyrF",
    "pyrK_2",
    "artM_1",
    "bioY2",
    "hmpT",
    "hypBA2",
    "EfmE745_03139",
    "EfmE745_03220",
    "EfmE745_03141",
    "EfmE745_03101",
    "EfmE745_03147",
    "EfmE745_03131",
    "EfmE745_01958",
    "EfmE745_01972",
    "EfmE745_02302",
    "EfmE745_00881",
    "rpsN2"
  ),

  ## “High-confidence, story-driving” genes (nice for clean summary plots)
  nucleotide_synthesis = c(
    # nucleotide
    "purD",
    "purH",
    "purL",
    "purQ",
    "purC",
    "purA",
    "guaB",
    "pyrF",
    "pyrK_2"
  ),
  pts = c(
    # PTS
    "manY_2",
    "manZ_3",
    "ptsL"
  ),
  cell_wall = c(
    # cell wall negative-fitness
    "clsA_1",
    "ddcP",
    "ldt_fmp",
    "mgs",
    "lytA_2"
  )
)

gff <- rtracklayer::import(file.path(
  dataFolder,
  grep(".gff", list.files(dataFolder), value = TRUE)
))
gff_df <- as.data.frame(mcols(gff)) |> as_tibble()

barplot_genes_fun <- function(
  genes_of_interest,
  genes_subset,
  subset_name,
  gff_df,
  logFC_df,
  heatmap_df,
  canu_sampleTable
) {
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

  ##################
  genes_of_interest <- genes_of_interest |>
    filter(Gene %in% genes_subset)

  logFC_genes_of_interest <- logFC_df_annotated |>
    filter(gene %in% genes_of_interest$Gene) |>
    left_join(genes_of_interest, by = c("gene" = "Gene"))
  # Cleanup duplicated columns and reorder
  logFC_genes_of_interest <- logFC_genes_of_interest |>
    # mutate(
    #   Log2FC = Log2FC.x,
    #   p_val = p_val.x,
    #   signif = signif.x
    # ) |>
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
  data_for_stats <- logFC_genes_of_interest |>
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
    mutate(condition = factor(condition, levels = c("serum", "bh")))

  # Keep only genes that have observations in both conditions
  valid_genes <- data_for_stats |>
    filter(!is.na(avg)) |>
    group_by(gene) |>
    summarize(n_conditions = n_distinct(condition), .groups = "drop") |>
    filter(n_conditions >= 2) |>
    pull(gene)

  if (length(valid_genes) > 0) {
    stats_df <- data_for_stats |>
      filter(gene %in% valid_genes) |>
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
  } else {
    stats_df <- tibble(
      gene = character(),
      group1 = character(),
      group2 = character(),
      y.position = numeric(),
      p = double(),
      p.adj = double(),
      p.signif = character(),
      p.adj.signif = character(),
      label = character(),
      stars = character()
    )
  }

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
    select(gene, Group, condition, sampleName, sample_count, avg, sd) |>
    mutate(
      avg_log10 = log10(avg),
      sample_count_log10 = log10(sample_count),
      ymin_log10 = log10(pmax(avg - sd, 1)),
      ymax_log10 = log10(avg + sd)
    )

  plot_y_values <- c(
    barplot_data$avg_log10,
    barplot_data$sample_count_log10,
    barplot_data$ymax_log10
  )
  plot_y_values <- plot_y_values[is.finite(plot_y_values)]
  if (length(plot_y_values) == 0) {
    plot_y_values <- c(0, 1)
  }
  y_limits <- range(plot_y_values)
  y_major_breaks <- seq(floor(y_limits[1]), ceiling(y_limits[2]), by = 1)
  if (length(y_major_breaks) == 0) {
    y_major_breaks <- 0:1
  }
  y_minor_breaks <- unlist(lapply(y_major_breaks, function(exp) {
    exp + log10(2:9)
  }))
  y_minor_breaks <- y_minor_breaks[
    y_minor_breaks > y_limits[1] & y_minor_breaks < y_limits[2]
  ]

  # Adust y.position for significance stars to be above the highest bar
  # Collapse barplot_data to one row per gene using the max(avg + sd)
  barplot_max <- barplot_data |>
    dplyr::group_by(gene) |>
    dplyr::summarise(
      y.position = {
        finite_y <- ymax_log10[is.finite(ymax_log10)]
        if (length(finite_y) == 0) NA_real_ else max(finite_y)
      },
      .groups = "drop"
    )

  # Join into stats_df
  stats_df <- stats_df |>
    select(-y.position) |>
    dplyr::left_join(barplot_max, by = "gene") |>
    mutate(
      y.position = dplyr::if_else(
        is.na(y.position),
        max(plot_y_values) + 0.12,
        y.position + 0.12
      )
    )

  # Compute actual p-values from normalized counts (fall back if rstatix returned NA)
  if (nrow(stats_df) > 0) {
    pvals_df <- normalized_counts_long |>
      group_by(gene) |>
      summarise(
        p_raw = tryCatch(
          t.test(count ~ condition)$p.value,
          error = function(e) {
            NA_real_
          }
        ),
        .groups = "drop"
      ) |>
      filter(gene %in% stats_df$gene)

    stats_df <- stats_df |>
      left_join(pvals_df, by = "gene") |>
      mutate(
        p = dplyr::coalesce(.data$p, .data$p_raw),
        p.signif = case_when(
          is.na(.data$p) ~ NA_character_,
          .data$p < 0.001 ~ "***",
          .data$p < 0.01 ~ "**",
          .data$p < 0.05 ~ "*",
          TRUE ~ NA_character_
        ),
        label = if_else(!is.na(p.signif), p.signif, stars)
      ) |>
      select(-p_raw)
  }

  # Build barplot
  plot_title <- genes_of_interest |>
    dplyr::distinct(Group) |>
    dplyr::pull(Group) |>
    na.omit() |>
    unique() |>
    paste(collapse = ", ")

  barplot_all <- ggplot(
    barplot_data,
    aes(x = gene, y = avg_log10, fill = condition)
  ) +
    geom_col(position = position_dodge(0.8), alpha = 0.8, width = 0.6) +
    geom_errorbar(
      aes(ymin = ymin_log10, ymax = ymax_log10),
      position = position_dodge(0.8),
      width = 0.2,
      linewidth = 0.5
    ) +
    geom_point(
      aes(y = sample_count_log10, color = condition),
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
    scale_y_continuous(
      labels = function(breaks) {
        scales::label_number()(10^breaks)
      },
      breaks = y_major_breaks,
      minor_breaks = y_minor_breaks
    ) +
    annotation_logticks(sides = "l", outside = FALSE) +
    labs(
      title = if (nzchar(plot_title)) plot_title else NULL,
      x = "Gene",
      y = "Normalized Counts",
      fill = "Condition"
    ) +
    # p-value annotations added conditionally below
    theme_classic() +
    theme(
      axis.ticks.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
      axis.text.y = element_text(size = 14),
      axis.title = element_text(size = 16),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14)
    )

  # Add p-value annotations before saving so they appear in the output image
  if (nrow(stats_df) > 0) {
    barplot_all <- barplot_all +
      ggpubr::stat_pvalue_manual(
        stats_df,
        label = "p.signif",
        y.position = "y.position",
        tip.length = 0
      )
  }

  ggsave(
    filename = paste0(
      "2026_VT\\Genome_analysis\\genome_analysis\\result\\barplot_",
      subset_name,
      ".png"
    ),
    plot = barplot_all,
    width = 14,
    height = 8,
    dpi = 300
  )
  return(barplot_all)
}

barplot_plots <- list()
for (i in seq_along(genes_subsets)) {
  barplot_plots[[i]] <- barplot_genes_fun(
    genes_of_interest = genes_of_interest,
    genes_subset = genes_subsets[[i]],
    subset_name = names(genes_subsets)[i],
    gff_df = gff_df,
    logFC_df = logFC_df,
    heatmap_df = heatmap_df,
    canu_sampleTable = canu_sampleTable
  )
}
