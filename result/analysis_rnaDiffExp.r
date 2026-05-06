# RNA differential expression analysis for the genome assembly project

# Libraries
library(DESeq2)
library(ggplot2)
library(ComplexHeatmap)
library(readr)

# Load data and assemble dataframe
dataFolder <- "2026_VT\\Genome_analysis\\genome_analysis\\result\\07_read_counts"

## From HTSeq-count (canu)
canu_files <- grep("canu", list.files(dataFolder), value = TRUE)

canu_condition <- sub("rna_(.*)_vs.*", "\\1", canu_files)

canu_sampleTable <- data.frame(
  sampleName = canu_condition,
  fileName = canu_files,
  condition = canu_condition
)

canu_sampleTable$condition <- factor(canu_sampleTable$condition)

ddsHTSeq_canu <- DESeqDataSetFromHTSeqCount(
  sampleTable = canu_sampleTable,
  directory = dataFolder,
  design = ~condition
)
# ddsHTSeq_canu

## From HTSeq-count (spades)
spades_files <- grep("spades", list.files(dataFolder), value = TRUE)
spades_condition <- sub("rna_(.*)_vs.*", "\\1", spades_files)
spades_sampleTable <- data.frame(
  sampleName = spades_condition,
  fileName = spades_files,
  condition = spades_condition
)
spades_sampleTable$condition <- factor(spades_sampleTable$condition)
ddsHTSeq_spades <- DESeqDataSetFromHTSeqCount(
  sampleTable = spades_sampleTable,
  directory = dataFolder,
  design = ~condition
)
# ddsHTSeq_spades

# Analysis & results
## canu
ddsHTSeq_canu_results <- DESeq(ddsHTSeq_canu)
res_canu <- results(ddsHTSeq_canu_results)
## spades
ddsHTSeq_spades_results <- DESeq(ddsHTSeq_spades)
res_spades <- results(ddsHTSeq_spades_results)
