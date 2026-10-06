# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

Workshop repository with two public bioinformatics datasets for hands-on analysis practice. No pre-written analysis scripts are included — the goal is to write analyses from scratch using real data.

## Data

- Raw data is in `data/raw/`. Never modify, move, or delete files there.
- Sample metadata: `data/raw/airway-rnaseq/samples.csv`

### airway-rnaseq

Bulk RNA-seq of human airway smooth muscle cells (Himes et al., PLoS One 2014; GEO GSE52778). 16 samples × 39,376 genes.

- `GSE52778_raw_counts_GRCh38.p13_NCBI.tsv.gz` — raw count matrix (NCBI Gene IDs)
- `Homo_sapiens.gene_info.gz` — gene ID → symbol/description mapping
- `samples.csv` — metadata: `sample_id`, `cell_line`, `treatment`, `dex`, `albuterol`

Design: 4 donors × 4 treatments (blocked design). Primary comparison: Untreated vs. Dexamethasone. ~9,300 genes are all-zero.

### moving-pictures-16s

16S rRNA amplicon sequencing from gut, palms, and tongue across two subjects over time (Caporaso et al., Genome Biology 2011; ENA ERP021896). 34 samples × 770 ASVs.

- `counts.tsv` — ASV count table
- `samples.tsv` — metadata: `body_site`, `subject`, `days_since_experiment_start`, `reported_antibiotic_usage`
- `taxonomy.tsv` — per-ASV taxonomy (Greengenes 13_8)

Sequencing depth varies widely (897–9,820 reads); the 5 shallowest samples are all right-palm from subject-1. Strongest signal is body site differences.

## Outputs

- Scripts go in `scripts/`, numbered in the order they run (e.g. `01_qc.R`, `02_deseq2.R`).
- All outputs go in `results/`. Figures go in `results/figures/` as PNG.
- Record software versions used in `results/versions.txt`.

## Style

- Explain statistical choices (normalization method, test, thresholds) in a comment at the top of each script.
