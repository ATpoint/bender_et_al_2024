# Bender et al 2024 EMBO

Code documentation for our paper **Bender et al (2024) EMBO**, 
[Redistribution of PU.1 partner transcription factor binding secures cell survival in leukemogenesis](https://www.embopress.org/doi/full/10.1038/s44318-024-00295-y).

## Content

- `00_packageVersions.Rmd`
Lists all R package and command line tool versions used in the Docker image we used for all analysis.

- `01_dataRetrieval.Rmd`
Uses biomaRt to fetch a lookup table of human-to-mouse orthologs necessary for some analysis scripts.

- `02_scRNAseq.Rmd`
All scRNA-seq analysis. Covers Figures 1A, 1B, 1C, 1D, 1E, 1F, 2A, 2B, 2C, Appendix Figure S1B/C/D/E/F/G/H/I/J and S4B.

- `03_RNAseq_Celllines.Rmd`
Bulk RNA-seq of cell lines. Covers Appendix Figure S2D.

- `04_shRNAscreen.Rmd`
The shRNA screen analysis. Covers Figures 3A, 3B, 3C and Appendix Figure S3B.

- `05_BeatAMLcohort.Rmd`
Integration of the BeatAML cohort. Covers Figure 4A and 4B.

- `06_proteome.Rmd`
Proteome analysis from Hox cell lines. Covers Figures 2D, 2E and 5A.

- `07_atacseq.Rmd`
ATAC-seq analysis and motif search. Covers Figures 6A-D, and 6I/J and data generation for manually drawn Figures 6E-H.

- `08_chipseq.Rmd`
RUNX1 ChIP-seq. Covers Figures 7A-E and Appendix Figure S5A/B.

- `09_other.Rmd`
The autohagy genes from RNA-seq in ex vivo murine cells and the low-throughput experiment analysis.
Covers Figure 3D, 4C, 5B-H, 7F/G, Appendix Figure S3F

## Run Code

Follow these steps to reproduce the figures:

1. Get source data

Download all preprocessed OMICS data (counts and metadata) from GEO, the numeric source data submitted to EMBO,
and the proteome dataset which we provide via this repository, since [the one submitted to EMBO](https://www.embopress.org/doi/suppl/10.1038/s44318-024-00295-y/suppl_file/44318_2024_295_moesm11_esm.xlsx) underwent some Excel-ish gene name to date conversion during their publication preparations.

```bash
mkdir -p source_data && cd source_data

# scRNA-seq
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE250nnn/GSE250629/suppl/GSE250629%5Fscrnaseq%5FrawCounts%5Funfiltered.mtx.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE250nnn/GSE250629/suppl/GSE250629%5Fscrnaseq%5Fcoldata%5Funfiltered.tsv.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE250nnn/GSE250629/suppl/GSE250629%5Fscrnaseq%5Frowdata%5Funfiltered.tsv.gz

# bulk RNA-seq from cell lines (don't be confused, B22 is the URE-AML cells, I forgot to relabel this before GEO submission)
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE250nnn/GSE250620/suppl/GSE250620%5Frnaseq%5Fcelllines%5FrawCounts.tsv.gz

# bulk RNA-seq from ex vivo cells
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE250nnn/GSE250622/suppl/GSE250622%5Frnaseq%5FexVivo%5FrawCounts.tsv.gz

# shRNA screen
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE250nnn/GSE250630/suppl/GSE250630%5Fshrna%5Fscreen%5FrawCounts.tsv.gz

# ATAC-seq
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE250nnn/GSE250619/suppl/GSE250619%5Fatacseq%5FrawCounts.tsv.gz

# ChIP-seq
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE251nnn/GSE251672/suppl/GSE251672%5Fchipseq%5Frunx1%5FrawCounts.tsv.gz

# Source data fo Figures 3,4,5,7 and the appendix figures
wget https://www.ebi.ac.uk/biostudies/files/S-SCDT-10_1038-S44318-024-00295-Y/Source_Data_Figure_3.zip
wget https://www.ebi.ac.uk/biostudies/files/S-SCDT-10_1038-S44318-024-00295-Y/Source_Data_Figure_4.zip
wget https://www.ebi.ac.uk/biostudies/files/S-SCDT-10_1038-S44318-024-00295-Y/Source_Data_Figure_5.zip
wget https://www.ebi.ac.uk/biostudies/files/S-SCDT-10_1038-S44318-024-00295-Y/Source_Data_Figure_7.zip
wget https://www.ebi.ac.uk/biostudies/files/S-SCDT-10_1038-S44318-024-00295-Y/Source_Data_Appendix.zip

# Unzip and only keep unzipped data
ls Source_Data_*.zip | while read p; do unzip ${p%.zip}; done
rm Source_Data_*.zip

# Proteome
wget https://github.com/ATpoint/bender_et_al_2024/raw/refs/heads/main/source_data/Dataset_EV_7.xlsx

# ChIP-seq peaks from published reanalyzed datasets. Code for how it was created it in the preprocessing documentation.
wget https://github.com/ATpoint/bender_et_al_2024/raw/refs/heads/main/source_data/LSK_PU1_IDR.txt.gz
wget https://github.com/ATpoint/bender_et_al_2024/raw/refs/heads/main/source_data/GMP_PU1_IDR.txt.gz
wget https://github.com/ATpoint/bender_et_al_2024/raw/refs/heads/main/source_data/GMP_CEBPA_IDR.txt.gz

# The tx2gene map that also contains the genomic coordinates of all TSS, made from the mouse GENCODE GTF file from version vM25 (Ensembl v100)
wget https://github.com/ATpoint/bender_et_al_2024/raw/refs/heads/main/source_data/tx2gene.txt.gz

# The subset of all relevant bigwig files for Figure 7D as RDS file, because the actual bigwigs are too big for easy sharing. 
# Email me if you need the bigwigs, I can provide.
wget https://github.com/ATpoint/bender_et_al_2024/raw/refs/heads/main/source_data/bigwig_signals.rds.xz
```

2. Run the Rmarkdown documents

We provide a Docker image that contains the exact software versions (R and command line applications) we used for analysis, allowing reproduction of the figures. Given that Docker is installed and in `PATH`, run this in your terminal:

Here, `DIR` is the full path to the directory with the folder `source_data` created above and the Rmarkdown and R scripts from the repository:

```r
DIR="/home/atpoint/bender_et_al"
IMAGE="atpoint/phd_project:1.9.5"
docker pull "$IMAGE" # takes some time, it's a big one due to legacy burden over many years...
docker run -d -p 8787:8787 -v "${DIR}":/projectdir -e PASSWORD=aVeryComplexPassword -e ROOT=TRUE -e IMAGE="$IMAGE" "$IMAGE"
```

Then type `localhost:8787` into your web browser with the username "rstudio" and password as above to access the interactive RStudio server session.
Be sure to render the code in chronological order since some scripts depend on output from scripts run before.

If not, or if you have any details questions, feel free to email me at `a.bender<guesswhat>uni-muenster.de` or open an issue.
