# Methods

## RNA-seq differential expression analysis

Raw read counts were obtained from NCBI GEO (accession GSE52778; NCBI reprocessing against GRCh38.p13), comprising 39,376 genes across 16 samples representing four primary human airway smooth muscle cell donors each treated with four conditions (Untreated, Dexamethasone, Albuterol, or Albuterol + Dexamethasone). Differential expression analysis was restricted to the eight samples comparing Untreated and Dexamethasone-treated cells (n = 4 donors per group).

Genes with fewer than 10 counts in fewer than four samples were excluded prior to modelling, retaining 16,802 genes for testing. Differential expression was performed using DESeq2 (v1.38.3) in R (v4.2.2) with a design formula of `~ cell_line + dex`, where `cell_line` (donor identity) was included as a blocking factor to account for the paired within-donor design, and `dex` (Untreated vs. Dexamethasone) was the variable of interest. Size factors were estimated by DESeq2's median-of-ratios method. Gene-wise dispersion estimates were fitted using a parametric mean-dispersion trend. Log2 fold changes were shrunk using the apeglm method (apeglm v1.20.0) to reduce noise for low-count genes without biasing estimates for highly expressed genes. The `|log2 fold change| > 1` significance threshold was applied to apeglm-shrunken values. Multiple testing correction used the Benjamini-Hochberg procedure; genes with a false discovery rate < 0.05 and |shrunken log2 fold change| > 1 were considered significant (453 upregulated, 420 downregulated in Dexamethasone relative to Untreated).

NCBI Gene IDs were converted to gene symbols using the NCBI *Homo sapiens* gene information table (Homo_sapiens.gene_info). For exploratory visualization, counts were variance-stabilized (VST, blind = TRUE for QC/PCA; blind = FALSE for the heatmap) using DESeq2. Principal component analysis was performed on the 500 most variable genes. Volcano plots were generated with ggplot2 (v4.0.2; Wickham 2016) and ggrepel (v0.9.6). The heatmap of the 30 most significant genes was produced with pheatmap (v1.0.13; Kolde 2019) on Z-scored VST values; columns were ordered Untreated then Dexamethasone with donor order held fixed within each group to preserve the paired structure. Additional R packages used: dplyr (v1.2.0), RColorBrewer (v1.1.3), scales (v1.4.0).

## Missing information

The following details a reviewer would likely ask about are not present in the analysis scripts or supporting files and could not be determined from the available data:

- **Read alignment software and version.** The count matrix was produced by NCBI's reprocessing pipeline; the specific aligner, version, and parameters used are not recorded in the local files.
- **Genome annotation version.** The genome build is stated as GRCh38.p13, but the gene annotation version (e.g., Ensembl release or RefSeq date) used during alignment/quantification is not recorded.
- **Read-level QC and trimming.** No trimming or adapter-removal step is present in the scripts; whether the NCBI reprocessing included these steps is not stated.
- **Sequencing platform and read length.** Not recorded in the local files (available from SRA but not pulled into the analysis).
- **ERCC spike-in handling.** The samples.csv file records which samples received ERCC spike-in mix 1; this information is not used in normalization or as a covariate, and no justification is given.
