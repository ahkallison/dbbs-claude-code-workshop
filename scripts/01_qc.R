# 01_qc.R — Quality control and exploratory analysis
#
# Design: 4 donors × 4 treatments (Untreated, Dexamethasone, Albuterol,
#   Albuterol+Dexamethasone) = 16 samples. Donor is a blocking factor.
#
# Normalization for visualization: Variance Stabilizing Transformation (VST)
#   via DESeq2. VST is preferred over raw/CPM for PCA because it homogenizes
#   variance across the mean-variance relationship inherent to count data.
#   VST is also faster and more numerically stable than rlog; DESeq2 recommends
#   rlog only for very small n (<12). With 16 samples here, VST is appropriate.
#
# Outputs:
#   results/figures/01a_library_sizes.png
#   results/figures/01b_genes_detected.png
#   results/figures/01c_pca.png
#   results/versions.txt  (written fresh; 02_deseq2.R will append to it)
#
# Required packages (install once):
#   install.packages(c("ggplot2","dplyr","RColorBrewer","scales"))
#   if (!requireNamespace("BiocManager")) install.packages("BiocManager")
#   BiocManager::install("DESeq2")

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(RColorBrewer)
  library(scales)
})

DATA  <- "data/raw/airway-rnaseq"
FIGS  <- "results/figures"
dir.create(FIGS, recursive = TRUE, showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# ── Publication theme ────────────────────────────────────────────────────────
pub_theme <- function() {
  theme_classic(base_size = 12, base_family = "Helvetica") +
    theme(
      plot.title    = element_text(face = "bold", size = 13, hjust = 0),
      plot.subtitle = element_text(size = 10, color = "grey40", hjust = 0),
      axis.text     = element_text(size = 10, color = "black"),
      axis.title    = element_text(size = 11),
      legend.title  = element_text(size = 10, face = "bold"),
      legend.text   = element_text(size = 9),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
      plot.margin   = margin(12, 16, 12, 12)
    )
}

# Palettes
TREAT_COLORS <- c(
  Untreated              = "#4393C3",
  Dexamethasone          = "#D6604D",
  Albuterol              = "#74C476",
  Albuterol_Dexamethasone = "#9970AB"
)
DONOR_SHAPES <- c(N61311 = 21, N052611 = 22, N080611 = 23, N061011 = 24)

# ── Load data ────────────────────────────────────────────────────────────────
message("Loading count matrix …")
counts <- read.table(
  gzfile(file.path(DATA, "GSE52778_raw_counts_GRCh38.p13_NCBI.tsv.gz")),
  header = TRUE, sep = "\t", row.names = 1, check.names = FALSE
)
meta <- read.csv(file.path(DATA, "samples.csv"), stringsAsFactors = FALSE)
meta$treatment <- factor(meta$treatment,
  levels = c("Untreated","Dexamethasone","Albuterol","Albuterol_Dexamethasone"))
meta$cell_line <- factor(meta$cell_line)

# Ensure column order matches metadata
counts <- counts[, meta$sample_id]
stopifnot(all(colnames(counts) == meta$sample_id))

message(sprintf("Matrix: %d genes × %d samples", nrow(counts), ncol(counts)))

# ── 1a. Library sizes ────────────────────────────────────────────────────────
lib_df <- data.frame(
  sample    = meta$sample_title,
  millions  = colSums(counts) / 1e6,
  treatment = meta$treatment,
  cell_line = meta$cell_line
)
lib_df$sample <- factor(lib_df$sample, levels = lib_df$sample)

p_lib <- ggplot(lib_df, aes(x = sample, y = millions, fill = treatment)) +
  geom_col(color = "white", linewidth = 0.3, width = 0.75) +
  geom_hline(yintercept = mean(lib_df$millions), linetype = "dashed",
             color = "grey40", linewidth = 0.5) +
  scale_fill_manual(values = TREAT_COLORS, name = "Treatment") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    title    = "Sequencing depth per sample",
    subtitle = sprintf("Dashed line = mean (%.1f M reads)", mean(lib_df$millions)),
    x = NULL, y = "Mapped reads (millions)"
  ) +
  pub_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

ggsave(file.path(FIGS, "01a_library_sizes.png"), p_lib,
       width = 9, height = 4.5, dpi = 300)
message("  Saved 01a_library_sizes.png")

# ── 1b. Genes detected (count > 0) ──────────────────────────────────────────
det_df <- data.frame(
  sample    = meta$sample_title,
  detected  = colSums(counts > 0),
  treatment = meta$treatment,
  cell_line = meta$cell_line
)
det_df$sample <- factor(det_df$sample, levels = det_df$sample)

p_det <- ggplot(det_df, aes(x = sample, y = detected, fill = treatment)) +
  geom_col(color = "white", linewidth = 0.3, width = 0.75) +
  geom_hline(yintercept = mean(det_df$detected), linetype = "dashed",
             color = "grey40", linewidth = 0.5) +
  scale_fill_manual(values = TREAT_COLORS, name = "Treatment") +
  scale_y_continuous(
    labels = scales::comma,
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    title    = "Genes detected per sample",
    subtitle = sprintf("Dashed line = mean (%s genes)", scales::comma(round(mean(det_df$detected)))),
    x = NULL, y = "Genes with ≥1 count"
  ) +
  pub_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

ggsave(file.path(FIGS, "01b_genes_detected.png"), p_det,
       width = 9, height = 4.5, dpi = 300)
message("  Saved 01b_genes_detected.png")

# ── 1c. PCA on VST-normalized counts ────────────────────────────────────────
message("Running VST and PCA …")
dds_all <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = meta,
  design    = ~ cell_line + treatment
)
# Filter low-count genes before VST.
# Threshold: ≥5 counts in ≥4 samples (out of 16 total). The DE script uses a
# stricter threshold (≥10 in ≥4) because it operates on the 8-sample subset
# only; with half the replicates, a higher count floor gives comparable
# sensitivity while reducing dispersion estimation noise.
# NOTE: samples.csv includes an ercc_mix column (ERCC spike-in mix 1 vs. none).
# N052611 and N061011 received mix 1; N61311 and N080611 received none. This
# batch variable is partially confounded with donor identity and is absorbed by
# the cell_line blocking term in the DE model. The PCA below should be examined
# to confirm it does not drive a visible batch axis independent of donor.
keep <- rowSums(counts(dds_all) >= 5) >= 4
dds_all <- dds_all[keep, ]
message(sprintf("  %d genes retained after filtering (count≥5 in ≥4 samples)", sum(keep)))

vst_all <- vst(dds_all, blind = TRUE)

# PCA on top 500 most-variable genes
pca_data <- plotPCA(vst_all, intgroup = c("treatment","cell_line"),
                    ntop = 500, returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"), 1)

pca_data$treatment <- factor(pca_data$treatment,
  levels = c("Untreated","Dexamethasone","Albuterol","Albuterol_Dexamethasone"))

p_pca <- ggplot(pca_data,
                aes(x = PC1, y = PC2, fill = treatment, shape = cell_line)) +
  geom_point(size = 4.5, stroke = 0.7, color = "black", alpha = 0.92) +
  scale_fill_manual(values = TREAT_COLORS, name = "Treatment",
                    labels = c("Untreated","Dexamethasone","Albuterol",
                               "Albuterol + Dex")) +
  scale_shape_manual(values = DONOR_SHAPES, name = "Donor") +
  guides(
    fill  = guide_legend(override.aes = list(shape = 21, size = 4)),
    shape = guide_legend(override.aes = list(fill = "grey60", size = 4))
  ) +
  labs(
    title    = "Principal component analysis of airway smooth muscle cells",
    subtitle = "VST-normalized counts, top 500 most-variable genes",
    x = sprintf("PC1  (%s%% variance)", pct_var[1]),
    y = sprintf("PC2  (%s%% variance)", pct_var[2])
  ) +
  pub_theme() +
  theme(legend.position = "right")

ggsave(file.path(FIGS, "01c_pca.png"), p_pca,
       width = 7.5, height = 5.5, dpi = 300)
message("  Saved 01c_pca.png")
message(sprintf("  PC1 = %.1f%%, PC2 = %.1f%%", pct_var[1], pct_var[2]))

# ── Versions ─────────────────────────────────────────────────────────────────
ver_lines <- c(
  paste0("# versions.txt — recorded ", Sys.time()),
  paste0("R: ", R.version$version.string),
  paste0("DESeq2: ", packageVersion("DESeq2")),
  paste0("ggplot2: ", packageVersion("ggplot2")),
  paste0("dplyr: ", packageVersion("dplyr")),
  paste0("RColorBrewer: ", packageVersion("RColorBrewer"))
)
writeLines(ver_lines, "results/versions.txt")

message("\nQC complete. Figures written to results/figures/")
