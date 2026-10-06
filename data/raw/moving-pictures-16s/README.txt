MOVING PICTURES 16S rRNA DATASET (QIIME 2 TUTORIAL SUBSET)

What this is
------------
16S rRNA gene amplicon sequencing of samples from two people, collected over
several months from four body sites: gut (stool), left palm, right palm and
tongue. Microbial communities differ strongly between body sites, so body
site is the main comparison.

This is the small subset of the original study that the QIIME 2 "Moving
Pictures" tutorial uses: 34 samples and 770 amplicon sequence variants
(ASVs, exact sequence types). The reads were processed with DADA2, a method
that removes sequencing errors to recover the exact sequences present
(single-end reads, trimmed to 120 bases). Each ASV was then classified using
the Greengenes 13_8 database, an older but widely used reference of
classified 16S sequences.

Study: Caporaso JG, Lauber CL, Costello EK, et al. Moving pictures of the
human microbiome. Genome Biology. 2011;12(5):R50.
DOI 10.1186/gb-2011-12-5-r50, PMID 21624126.

Raw data: ENA study ERP021896 (BioProject PRJEB19825), Qiita study 550.


Files
-----
counts.tsv
    Read counts, 770 ASVs (rows) x 34 samples (columns). Tab-separated.
    The first column, feature_id, is the ASV identifier (an MD5 hash of the
    ASV's DNA sequence). The other columns are sample IDs. Values are
    whole-number read counts.

taxonomy.tsv
    One row per ASV with its classification: feature_id, kingdom, phylum,
    class, order, family, genus, species, confidence. An empty cell means the
    classifier could not assign that rank with confidence. Names follow
    Greengenes 13_8, so some are older names (for example Bacteroidetes and
    Firmicutes rather than Bacteroidota and Bacillota), and names in square
    brackets, like [Prevotella], are provisional Greengenes names.

samples.tsv
    One row per sample.
      sample_id                    matches the count column names
      barcode_sequence             sequencing barcode
      body_site                    gut, left palm, right palm, tongue
      year, month, day             collection date
      subject                      subject-1 or subject-2
      reported_antibiotic_usage    Yes or No
      days_since_experiment_start  0, 84, 112, 140 or 168

    Samples per body site: gut 8, left palm 8, right palm 9, tongue 9.

original/
    The three files as published by the QIIME 2 tutorial, unmodified:
    table.qza, taxonomy.qza and sample-metadata.tsv. A .qza file is a zip
    archive. Inside, the count table is a BIOM file (an HDF5 binary format)
    and the taxonomy is a TSV. counts.tsv, taxonomy.tsv and samples.tsv were
    made from these files.


Where the files came from
-------------------------
  https://moving-pictures-tutorial.readthedocs.io/en/stable/data/moving-pictures/table.qza
  https://moving-pictures-tutorial.readthedocs.io/en/stable/data/moving-pictures/taxonomy.qza
  https://moving-pictures-tutorial.readthedocs.io/en/stable/data/moving-pictures/sample-metadata.tsv


Things to know before you analyze
---------------------------------
- Sequencing depth varies from 897 to 9,820 reads per sample. The five
  shallowest samples are all right palm samples from subject-1. Account for
  depth before comparing diversity between samples.
- Body site is the strongest comparison. Subject is the second.
- Look at how reported_antibiotic_usage relates to days_since_experiment_start
  before you test for an antibiotic effect.
