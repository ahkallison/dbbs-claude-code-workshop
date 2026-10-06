# 02_deseq2.R — Differential expression: Untreated vs. Dexamethasone
#
# Subset: 8 samples (4 donors × Untreated/Dexamethasone only).
#
# Method: DESeq2 (Love et al., Genome Biology 2014).
#   Chosen because:
#   (1) Counts are overdispersed integers — the negative binomial model fitted
#       by DESeq2 is better calibrated than normal/Poisson models.
#   (2) DESeq2's median-of-ratios size factors correct for differences in
#       sequencing depth without the distributional assumptions of TPM/RPKM.
#   (3) apeglm LFC shrinkage (Zhu et al. 2019) reduces noise for low-count
#       genes without biasing high-count genes, giving a better-ranked list.
#
# Design: ~ cell_line + dex
#   Donor (cell_line) is a blocking factor — this is a within-donor comparison.
#   'dex' encodes the treatment (no = Untreated, yes = Dexamethasone).
#
# Multiple testing: Benjamini-Hochberg FDR (DESeq2 default).
# Significance thresholds: padj < 0.05, |log2FC| > 1.
#
# Gene ID conversion: NCBI Homo_sapiens.gene_info (GeneID → Symbol).
#
# Outputs:
#   results/dex_vs_untreated_all.csv
#   results/analysis_summary.txt
#   results/figures/02a_volcano.png
#   results/figures/02b_heatmap.png
#   results/versions.txt  (appended to file written by 01_qc.R)
#
# Required packages (install once):
#   install.packages(c("ggplot2","ggrepel","dplyr","RColorBrewer","scales"))
#   if (!requireNamespace("BiocManager")) install.packages("BiocManager")
#   BiocManager::install(c("DESeq2","apeglm","pheatmap"))

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(dplyr)
  library(RColorBrewer)
  library(scales)
})

DATA <- "data/raw/airway-rnaseq"
FIGS <- "results/figures"
dir.create(FIGS, recursive = TRUE, showWarnings = FALSE)

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

# ── Known glucocorticoid-response genes ──────────────────────────────────────
GC_GENES <- c("TSC22D3","PER1","DUSP1","KLF15","FKBP5","CRISPLD2")

# ── Load data ────────────────────────────────────────────────────────────────
message("Loading counts and metadata …")
counts <- read.table(
  gzfile(file.path(DATA, "GSE52778_raw_counts_GRCh38.p13_NCBI.tsv.gz")),
  header = TRUE, sep = "\t", row.names = 1, check.names = FALSE
)
meta <- read.csv(file.path(DATA, "samples.csv"), stringsAsFactors = FALSE)

# Subset to Untreated and Dexamethasone only
meta_sub <- meta[meta$treatment %in% c("Untreated","Dexamethasone"), ]
counts_sub <- counts[, meta_sub$sample_id]
meta_sub$dex      <- factor(meta_sub$dex, levels = c("no","yes"))
meta_sub$cell_line <- factor(meta_sub$cell_line)

message(sprintf("Subset: %d samples, %d genes", ncol(counts_sub), nrow(counts_sub)))

# ── Gene symbol mapping ───────────────────────────────────────────────────────
message("Loading gene symbol table …")
gene_info <- read.table(
  gzfile(file.path(DATA, "Homo_sapiens.gene_info.gz")),
  header = TRUE, sep = "\t", quote = "", comment.char = "",
  fill = TRUE, stringsAsFactors = FALSE
)
# Keep only columns we need; column names vary slightly across releases
id_col  <- grep("GeneID|Gene_ID",  names(gene_info), value = TRUE)[1]
sym_col <- grep("^Symbol$",        names(gene_info), value = TRUE)[1]
desc_col<- grep("description|full_name", names(gene_info), ignore.case=TRUE, value=TRUE)[1]

gene_map <- gene_info[, c(id_col, sym_col, desc_col)]
names(gene_map) <- c("GeneID","Symbol","Description")
gene_map <- gene_map[!duplicated(gene_map$GeneID), ]
rownames(gene_map) <- as.character(gene_map$GeneID)

message(sprintf("  %d GeneIDs in map", nrow(gene_map)))

# ── DESeq2 ───────────────────────────────────────────────────────────────────
message("Running DESeq2 …")
dds <- DESeqDataSetFromMatrix(
  countData = counts_sub,
  colData   = meta_sub,
  design    = ~ cell_line + dex
)

# Pre-filter: keep genes with ≥10 counts in at least 4 samples (= all replicates
# of one group). Threshold is stricter than the QC script (which uses ≥5 on all
# 16 samples) because the DE subset has only 8 samples; a higher count floor
# gives comparable sensitivity while reducing dispersion estimation noise.
# NOTE: samples.csv includes an ercc_mix column — donors N052611 and N061011
# received ERCC mix 1; N61311 and N080611 received none. Because ercc_mix is
# partially confounded with cell_line, the blocking term in the design absorbs
# most of this batch effect, but it is not explicitly modelled.
keep <- rowSums(counts(dds) >= 10) >= 4
dds  <- dds[keep, ]
message(sprintf("  %d genes after pre-filtering", sum(keep)))

dds  <- DESeq(dds, quiet = TRUE)

# LFC shrinkage with apeglm (better than normal/ashr for ranking)
res_shrunk <- lfcShrink(dds, coef = "dex_yes_vs_no", type = "apeglm", quiet = TRUE)
res_df <- as.data.frame(res_shrunk)
res_df$GeneID <- rownames(res_df)

# Add raw (unshrunken) LFC for reference
res_raw <- results(dds, contrast = c("dex","yes","no"))
res_df$log2FC_raw <- res_raw[rownames(res_df), "log2FoldChange"]

# Merge gene symbols
res_df$Symbol <- gene_map[res_df$GeneID, "Symbol"]
res_df$Description <- gene_map[res_df$GeneID, "Description"]

# Replace NA symbols with GeneID
res_df$Symbol[is.na(res_df$Symbol)] <- res_df$GeneID[is.na(res_df$Symbol)]

# Sort by adjusted p-value
res_df <- res_df[order(res_df$padj, na.last = TRUE), ]

# Write full results
write.csv(res_df, "results/dex_vs_untreated_all.csv", row.names = FALSE)
message(sprintf("  Results written: %d genes", nrow(res_df)))

# Summary. n per group recorded here for provenance.
n_untreated   <- sum(meta_sub$dex == "no")
n_dex         <- sum(meta_sub$dex == "yes")
message(sprintf("  Group sizes: Untreated n=%d, Dexamethasone n=%d", n_untreated, n_dex))

# |LFC| > 1 threshold is applied to the apeglm-shrunken log2FoldChange column
# (res_df$log2FoldChange), not to log2FC_raw (the unshrunken MLE). Shrunken
# values are used throughout for consistency with the volcano plot.
sig <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, ]
message(sprintf("  Significant (padj<0.05, |shrunken LFC|>1): %d up, %d down",
                sum(sig$log2FoldChange > 0), sum(sig$log2FoldChange < 0)))

# ── Known gene check ─────────────────────────────────────────────────────────
message("\nKnown glucocorticoid-response genes:")
known <- res_df[res_df$Symbol %in% GC_GENES, c("Symbol","log2FoldChange","padj")]
known <- known[order(known$padj, na.last = TRUE), ]
print(known, row.names = FALSE)

# ── Volcano plot ─────────────────────────────────────────────────────────────
message("\nGenerating volcano plot …")

vol_df <- res_df
vol_df <- vol_df[!is.na(vol_df$padj) & !is.na(vol_df$log2FoldChange), ]

# Significance categories
vol_df$category <- "NS"
vol_df$category[vol_df$padj < 0.05 & vol_df$log2FoldChange >  1] <- "Up"
vol_df$category[vol_df$padj < 0.05 & vol_df$log2FoldChange < -1] <- "Down"
vol_df$category <- factor(vol_df$category, levels = c("Up","Down","NS"))

# Cap -log10(padj) for display (some padj are effectively 0)
vol_df$neg_log10p <- pmin(-log10(vol_df$padj), 320)

# Flag known genes
vol_df$is_known <- vol_df$Symbol %in% GC_GENES
label_df <- vol_df[vol_df$is_known, ]

cat_colors <- c(Up = "#D6604D", Down = "#4393C3", NS = "grey75")
cat_alpha  <- c(Up = 0.75, Down = 0.75, NS = 0.35)
cat_size   <- c(Up = 1.2,  Down = 1.2,  NS = 0.6)

p_vol <- ggplot(vol_df, aes(x = log2FoldChange, y = neg_log10p)) +
  # All points
  geom_point(aes(color = category, alpha = category, size = category)) +
  # Known gene overlay
  geom_point(data = label_df, aes(fill = category),
             shape = 21, color = "black", size = 3, stroke = 0.8) +
  geom_label_repel(
    data          = label_df,
    aes(label = Symbol),
    size          = 3.2,
    fontface      = "bold.italic",
    box.padding   = 0.5,
    point.padding = 0.4,
    segment.color = "grey40",
    segment.size  = 0.4,
    fill          = "white",
    alpha         = 0.9,
    max.overlaps  = 20
  ) +
  # Threshold lines
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  scale_color_manual(values = cat_colors, name = NULL,
    labels = c(
      Up   = sprintf("Up in Dex (%d)", sum(vol_df$category=="Up")),
      Down = sprintf("Down in Dex (%d)", sum(vol_df$category=="Down")),
      NS   = sprintf("NS (%d)", sum(vol_df$category=="NS"))
    )) +
  scale_fill_manual(values = cat_colors, guide = "none") +
  scale_alpha_manual(values = cat_alpha, guide = "none") +
  scale_size_manual(values = cat_size, guide = "none") +
  scale_x_continuous(limits = c(-8.5, 8.5), breaks = seq(-8, 8, 2)) +
  labs(
    title    = "Dexamethasone vs. Untreated — airway smooth muscle cells",
    subtitle = "DESeq2 with apeglm LFC shrinkage  |  Design: ~ donor + dexamethasone  |  4 donors, paired",
    x        = expression(log[2]~"fold change (Dex / Untreated)"),
    y        = expression(-log[10]~"(adjusted"~italic(p)*"-value)")
  ) +
  annotate("text", x = 5.5, y = -log10(0.05) + 6, label = "FDR = 5%",
           size = 3, color = "grey40") +
  pub_theme() +
  theme(legend.position = c(0.85, 0.90),
        legend.background = element_rect(fill = alpha("white", 0.8), color = NA))

ggsave(file.path(FIGS, "02a_volcano.png"), p_vol,
       width = 8, height = 6.5, dpi = 300)
message("  Saved 02a_volcano.png")

# ── Heatmap of top 30 genes ──────────────────────────────────────────────────
message("Generating heatmap …")

# VST on the 8-sample subset for heatmap
dds_vst <- estimateSizeFactors(dds)
vst_sub  <- vst(dds_vst, blind = FALSE)
vst_mat  <- assay(vst_sub)

# Top 30 by adjusted p-value (among genes with symbol mapping)
top30_ids <- head(res_df$GeneID[!is.na(res_df$padj) &
                                 !is.na(res_df$Symbol) &
                                 res_df$Symbol != res_df$GeneID], 30)
top30_mat  <- vst_mat[top30_ids, ]
top30_syms <- res_df$Symbol[match(top30_ids, res_df$GeneID)]
rownames(top30_mat) <- top30_syms

# Z-score rows
top30_scaled <- t(scale(t(top30_mat)))

# Order columns: Untreated then Dexamethasone, donor order fixed within each
# group. Column clustering is disabled so the paired structure that the DE model
# uses is preserved visually — free clustering would obscure within-donor pairing.
col_order <- order(meta_sub$dex, meta_sub$cell_line)
top30_scaled <- top30_scaled[, col_order]
meta_ordered <- meta_sub[col_order, ]

col_ann <- data.frame(
  Treatment = ifelse(meta_ordered$dex == "yes", "Dexamethasone", "Untreated"),
  Donor     = meta_ordered$cell_line,
  row.names = meta_ordered$sample_id
)
colnames(top30_scaled) <- meta_ordered$sample_title

ann_colors <- list(
  Treatment = c(Dexamethasone = "#D6604D", Untreated = "#4393C3"),
  Donor     = c(N61311 = "#F6A800", N052611 = "#3E8E41",
                N080611 = "#7B3FA0", N061011 = "#1A6FAF")
)

png(file.path(FIGS, "02b_heatmap.png"),
    width = 8, height = 9, units = "in", res = 300)
pheatmap(
  top30_scaled,
  color            = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  breaks           = seq(-3, 3, length.out = 101),
  cluster_rows     = TRUE,
  cluster_cols     = FALSE,
  annotation_col   = col_ann,
  annotation_colors = ann_colors,
  show_colnames    = TRUE,
  show_rownames    = TRUE,
  fontsize_row     = 9,
  fontsize_col     = 8,
  fontsize         = 9,
  border_color     = NA,
  main             = "Top 30 DEGs: Dexamethasone vs. Untreated\n(VST-normalized, Z-scored per gene; columns: Untreated > Dex, paired by donor)"
)
dev.off()
message("  Saved 02b_heatmap.png")

# ── Write complete versions file (overwrites to avoid duplicates on reruns) ───
ver_lines <- c(
  paste0("# versions.txt — recorded ", Sys.time()),
  paste0("R: ", R.version$version.string),
  paste0("DESeq2: ",       packageVersion("DESeq2")),
  paste0("apeglm: ",       packageVersion("apeglm")),
  paste0("ggplot2: ",      packageVersion("ggplot2")),
  paste0("ggrepel: ",      packageVersion("ggrepel")),
  paste0("pheatmap: ",     packageVersion("pheatmap")),
  paste0("dplyr: ",        packageVersion("dplyr")),
  paste0("RColorBrewer: ", packageVersion("RColorBrewer")),
  paste0("scales: ",       packageVersion("scales"))
)
writeLines(ver_lines, "results/versions.txt")

# Record group sizes in a self-contained summary file
summary_lines <- c(
  "# analysis_summary.txt",
  paste0("comparison: Dexamethasone vs. Untreated"),
  paste0("n_untreated: ", n_untreated),
  paste0("n_dexamethasone: ", n_dex),
  paste0("n_donors: ", nlevels(meta_sub$cell_line)),
  paste0("design: paired (blocked by donor)"),
  paste0("genes_tested: ", nrow(res_df)),
  paste0("significant_padj05_lfc1_up: ", sum(sig$log2FoldChange > 0)),
  paste0("significant_padj05_lfc1_down: ", sum(sig$log2FoldChange < 0)),
  paste0("lfc_column_used_for_threshold: apeglm-shrunken (log2FoldChange)"),
  paste0("generated: ", Sys.time())
)
writeLines(summary_lines, "results/analysis_summary.txt")

message("\nDESeq2 analysis complete.")
