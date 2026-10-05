AIRWAY RNA-SEQ DATASET (GEO GSE52778)

What this is
------------
Bulk RNA-seq of primary human airway smooth muscle cells from 4 donors. Each
donor's cells were left untreated or treated with dexamethasone (a
glucocorticoid steroid), albuterol (a bronchodilator), or both, for 16 samples
in total. Dexamethasone changes the expression of hundreds of genes in these
cells, so the untreated vs dexamethasone comparison gives a strong, easy-to-see
signal.

Study: Himes BE, Jiang X, Wagner P, et al. RNA-Seq transcriptome profiling
identifies CRISPLD2 as a glucocorticoid responsive gene that modulates cytokine
function in airway smooth muscle cells. PLoS One. 2014;9(6):e99625.
DOI 10.1371/journal.pone.0099625, PMID 24926665.

GEO series: GSE52778 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778)
SRA study: SRP033351. Runs SRR1039508 to SRR1039523, one per sample.


Files
-----
GSE52778_raw_counts_GRCh38.p13_NCBI.tsv.gz
    Raw read counts per gene, 39,376 genes x 16 samples. Tab-separated.
    The first column, GeneID, is the NCBI Gene ID (a number, for example 1831
    is TSC22D3). The other columns are GEO sample IDs (GSM...). Values are
    whole-number counts. NCBI produced this table by reprocessing the raw
    reads against the GRCh38.p13 human genome; it is not the authors' own
    count table.

samples.csv
    One row per sample, in the same order as the count columns.
      sample_id     GEO sample ID, matches the count column names
      sample_title  the label the authors gave the sample
      cell_line     donor cell line (4 donors)
      treatment     Untreated, Dexamethasone, Albuterol, or
                    Albuterol_Dexamethasone
      dex           yes/no, whether the treatment included dexamethasone
      albuterol     yes/no, whether the treatment included albuterol
      ercc_mix      ERCC spike-in mix used, if any ("-" means none)
      srr, srx, biosample   SRA and BioSample IDs for the raw reads
    Built from the GEO series matrix file and SRA run information.

GSE52778_series_matrix.txt.gz
    GEO's original sample description file, unmodified. samples.csv was
    built from it.

Homo_sapiens.gene_info.gz
    NCBI's human gene table, for turning GeneIDs into gene symbols and
    descriptions. Tab-separated; the columns you need are GeneID, Symbol and
    description. NCBI updates this file daily, and about 1,700 of the GeneIDs
    in the count table are missing from the current release.


Where the files came from
-------------------------
If you want to download them yourself:

  Counts:
  https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=GSE52778&format=file&file=GSE52778_raw_counts_GRCh38.p13_NCBI.tsv.gz

  Series matrix:
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE52nnn/GSE52778/matrix/GSE52778_series_matrix.txt.gz

  Gene table:
  https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz

GEO's download page sometimes answers scripted downloads with a "checking
your browser" page instead of the file. If a downloaded .gz file will not
open, check it with: file <name>.gz
It should say "gzip compressed data".


Things to know before you analyze
---------------------------------
- Each donor received every treatment, so donor is a blocking factor:
  compare treatments within donor, not across donors.
- A good first comparison is Untreated vs Dexamethasone (8 samples, 4 per
  group).
- Library sizes range from about 19 to 39 million counted reads, so
  normalize before comparing samples.
- About 9,300 genes have zero counts in every sample.
