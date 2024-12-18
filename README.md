! code will be available by tomorrow (19th december) -- sorry for the delay !

# Bender et al 2024 EMBO

This repository contains the code documentation for our paper **Bender et al (2024) EMBO**, 
[Redistribution of PU.1 partner transcription factor binding secures cell survival in leukemogenesis](https://www.embopress.org/doi/full/10.1038/s44318-024-00295-y).

It reproduces all figures based on either the GEO-submitted preprocessed OMICS data or the EMBO-submitted numerical source data.

In case you want to hands-on run the analysis, follow below steps:

## Get Source Data Data

Download all preprocessed OMICS data (counts and metadata) from GEO:

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
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE250nnn/GSE250630/suppl/GSE250630%5Fshrna%5Fscreen%5FrawCounts.tsv.gz

# ATAC-seq
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE250nnn/GSE250619/suppl/GSE250619%5Fatacseq%5FrawCounts.tsv.gz

# ChIP-seq
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE251nnn/GSE251672/suppl/GSE251672%5Fchipseq%5Frunx1%5FrawCounts.tsv.gz

# Source data fo Figures 3,4,5,7 and the appendix figures
wget https://www.ebi.ac.uk/biostudies/files/S-SCDT-10_1038-S44318-024-00295-Y/Source_Data_Figure_3.zip
wget https://www.ebi.ac.uk/biostudies/files/S-SCDT-10_1038-S44318-024-00295-Y/Source_Data_Figure_4.zip
wget https://www.ebi.ac.uk/biostudies/files/S-SCDT-10_1038-S44318-024-00295-Y/Source_Data_Figure_5.zip
wget https://www.ebi.ac.uk/biostudies/files/S-SCDT-10_1038-S44318-024-00295-Y/Source_Data_Figure_7.zip
wget https://www.ebi.ac.uk/biostudies/files/S-SCDT-10_1038-S44318-024-00295-Y/Source_Data_Appendix.zip

# Unzip and only keep unzipped data
ls Source_Data_*.zip | while read p; do unzip ${p%.zip}; done
rm Source_Data_*.zip
```

The raw (that is MaxQuant-processed) proteome is available as an Excel sheet via the EMBO website ([here -- don't use](https://www.embopress.org/doi/suppl/10.1038/s44318-024-00295-y/suppl_file/44318_2024_295_moesm11_esm.xlsx)) as an Appendix file, but the some gene names were converted to dates during the submission process, so we provide the "proper" spreadsheet via
this GitHub repository:

```bash
wget https://github.com/ATpoint/bender_et_al_2024/raw/refs/heads/main/Dataset_EV_7.xlsx
```

## Run Full Analysis

We provide a Docker image that contains the exact software versions (R and command line applications) we used for analysis, allowing reproduction of the figures. Given that Docker is installed and in `PATH`, run this in your terminal:

Here, `DIR` is the full path to the directory with the downloaded source data from above. 

```r
DIR="/path/to/bender_et_al_2024/"
IMAGE="atpoint/phd_project:1.9.5"
docker pull "$IMAGE" # takes some time, it's a big one due to legacy burden over many years...
docker run -d -p 8787:8787 -v "${DIR}":/projectdir -e PASSWORD=aVeryComplexPassword -e ROOT=TRUE -e IMAGE="$IMAGE" "$IMAGE"
```

Then type `localhost:8787` into your web browser with the username "rstudio" and password as above to access the interactive RStudio server session. If you run on Codespaces you can forward the 8787 port to your local machine and access the session that way.

From here, run the Rmarkdown scripts in chronological order.

If not, or if you have any details questions, feel free to email me at `a.bender<guesswhat>uni-muenster.de` or open an issue.
