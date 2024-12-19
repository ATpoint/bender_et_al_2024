# Run before execution of the Rmarkdown documents

knitr::opts_chunk$set(
  eval = TRUE,
  message = TRUE,
  warning = TRUE,
  comment = "#>"
)

# Define rootdir and cachedir
outdir <- paste0(rootdir, "/analysis_results/")
suppressWarnings(dir.create(outdir, recursive = TRUE))

# Disable BLAS multithreading
RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)

# For MAGeCK
Sys.setenv(MAGECKPATH = "/home/rstudio/.local/bin/")

suppressMessages({
  library(AnnotationDbi)
  library(ashr)
  library(batchelor)
  library(BiocParallel)
  library(biomaRt)
  library(Biostrings)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(chromVAR)
  library(ComplexHeatmap)
  ht_opt$message <- FALSE
  library(CreateGeneSignatures)
  library(data.table)
  library(DESeq2)
  library(edgeR)
  library(EnrichedHeatmap)
  library(fgsea)
  library(GenomicRanges)
  library(GEOquery)
  library(ggpointdensity)
  library(ggpubr)
  library(ggrepel)
  library(ggsignif)
  library(ggupset)
  library(glmGamPoi)
  library(gprofiler2)
  library(gridtext)
  library(hexbin)
  library(imputeLCMD)
  library(limma)
  library(magrittr)
  library(matrixStats)
  library(Matrix)
  library(memes)
  library(motifmatchr)
  library(MSnbase)
  library(openxlsx)
  library(org.Mm.eg.db)
  library(patchwork)
  library(S4Vectors)
  library(scater)
  library(scDblFinder)
  library(scran)
  library(sctransform)
  library(scuttle)
  library(SingleCellExperiment)
  library(SingleR)
  library(slingshot)
  library(SummarizedExperiment)
  library(TFBSTools)
  library(tidyverse)
  library(tximport)
  library(viridis)
  library(vizzy)
  library(UCell)
})

source(paste0(rootdir, "/functions.R"))

# BiocParallel (always use bplapply instead of mclapply)
mc_workers <- floor(parallel::detectCores() * 0.9)
doParallel::registerDoParallel(mc_workers)
bpparam <- BiocParallel::DoparParam()

list.ggplot <- list()

list.ggplot$colorblind_cols <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#661100",
  "#332288", "gray40", "black"
)

# make brighter and darker
list.ggplot$colorblind_cols <-
  c(
    list.ggplot$colorblind_cols,
    colorspace::lighten(list.ggplot$colorblind_cols, .5),
    colorspace::lighten(list.ggplot$colorblind_cols, .9)
  )

# some global options for ggplot, to achieve some standard appearances:
list.ggplot$theme <- ggplot2::theme_bw # default theme
list.ggplot$themesize <- 12 # base_size option
list.ggplot$pointsize <- 0.5 # point size for UMAPs
list.ggplot$pointsize_larger <- 3 # for dotplots with few points
list.ggplot$textsize <- 3.5 #
list.ggplot$ablinesize <- .5 #
list.ggplot$linesize <- 1 #
list.ggplot$umap_alpha <- 0.6 # for UMAPs: alpha of points
list.ggplot$legendsize <- 3
list.ggplot$genotype_colors <- c("#999933", "#4d77f5")

# set the global theme:
theme_set(list.ggplot$theme(base_size = list.ggplot$themesize))
theme_update(
  axis.text = element_text(size = list.ggplot$themesize),
  axis.title = element_text(size = list.ggplot$themesize),
  legend.text = element_text(size = list.ggplot$themesize),
  legend.title = element_text(size = list.ggplot$themesize),
  strip.text = element_text(size = list.ggplot$themesize),
  strip.background = element_rect(fill = "white"),
  text = element_text(size = list.ggplot$themesize)
)

gg.noX <- theme(
  axis.title.x = element_blank(),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank()
)

base::options(
  ggplot2.discrete.colour = list.ggplot$colorblind_cols,
  ggplot2.discrete.fill = list.ggplot$colorblind_cols
)