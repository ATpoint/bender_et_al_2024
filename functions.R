####################
# CUSTOM FUNCTIONS #
####################

#------------------------------------
# Percent of a group expressing each gene
#------------------------------------
get_pexpr <- function(data, group, threshold=0, digits=2){
  
  if(ncol(data)!=length(group)) stop("ncol(data) != length(group)")
  if(!is.numeric(threshold) | threshold < 0) stop("threshold must be numeric and > 0")
  
  datar <- (data>threshold) * 1
  a <- base::rowsum(x=t(datar), group=group)
  b <- as.numeric(table(group)[rownames(a)])
  f <- base::round(100*t(apply(a, 2, function(x) x/b)), digits=digits)
  f
  
}

#------------------------------------
# Plotting UMAPs
#------------------------------------

#' Plot UMAPs with ggplot
#' 
#' A function to plot UMAPs based on a SingleCellExperiment accessing the 
#' coordinates in the reducedDims() slot coloring either by colData features
#' or by gene expression
#' 
#' @param sce a SingleCellExperiment
#' @param gene a gene name to extract expr values for, assumed to be "genename" when rownames(sce) are geneid_genename,
#' @param gene_col the color for points with high expression values, low expr is hardcoded as grey
#' @param by use this column from colData for coloring
#' @param text_add logical, whether to add the by label to the plot, e.g. the name of the cluster
#' @param text_col color for text_add
#' @param text_use_label logical, whether to add text as label box rather than just letters (geom_label instead of geom_text)
#' @param legend_title guess what
#' @param embed_onto name of the reducedDim to use as manifold
#' @param embed_dims vector with two integers, the column idx of embed_onto, e.g. c(1,2) when using the first and second UMAP dimensions
#' @param return_data logical, whether to return a list with the plot and the data.frame with the underlying data
plot_umap <- function(sce, 
                      gene=NULL, 
                      gene_col="darkblue",
                      by=NULL, text_add=TRUE, text_col="black", text_use_label=FALSE,
                      legend_title=NULL,
                      embed_onto="UMAP", embed_dims=c(1,2), # embed_dims are the columns or reducedDims, e.g. UMAP1+2
                      return_data=FALSE, label_size=NULL){
  
  suppressMessages({
    require(SingleCellExperiment)
    require(scater)
    require(tidyverse)
    require(viridis)
  })
  
  if(!is(sce, "SingleCellExperiment"))
    stop("sce must be a SingleCellExperiment")
  
  # check that only by or gene is used:
  if(!is.null(gene) & !is.null(by))
    stop("Use either by or gene, not both!")
  
  # the dim reduction to use as backbone for the plot:
  if(!embed_onto %in% names(reducedDims(sce)))
    stop("embed_onto is not a reducedDim in the sce!")
  
  colData(sce) <- droplevels.data.frame(colData(sce))
  
  plotobj <- as.data.frame(reducedDim(sce, embed_onto))
  #rm(sce) # otherwise this will be part of the ggplot environment!
  
  if(!all(embed_dims %in% 1:ncol(plotobj)))
    stop("Check embed_dims -- some are out of bounds")
  plotobj <- plotobj[,embed_dims]
  colnames(plotobj) <- paste0("dim", embed_dims)
  
  # plot gene expression:
  if(!is.null(gene)){
    
    grepped <- grep(paste0("_", gene, "$"), rownames(sce), value=TRUE)
    if(length(grepped)==0)
      stop("gene not found!")
    
    plotobj$values <- as.numeric(retrieveCellInfo(sce, grepped)$value)
    
    if(is.null(legend_title)) legend_title <- gene
    
    ggobj <- 
      plotobj %>%
      arrange(values) %>%
      ggplot(aes(x=dim1, y=dim2, color=values)) +
      geom_point(alpha=list.ggplot$umap_alpha, size=list.ggplot$pointsize) +
      scale_color_viridis(name=gene)
    
  }
  
  # plot colData elements:
  if(!is.null(by)){
    
    retrieved_by <- tryCatch(expr=retrieveCellInfo(sce, by),
                             error=function(e) stop("by not found!"))
    
    if(class(retrieved_by$value) != "numeric"){
      
      plotobj$group <- factor(retrieved_by$value)  
      
      text_df <- 
        sapply(levels(plotobj$group), function(x){
          data.frame(group=x,
                     x=median(plotobj[plotobj$group==x,"dim1"]),
                     y=median(plotobj[plotobj$group==x,"dim2"]))
        }, simplify=FALSE) 
      
      text_df <- do.call(rbind, text_df)
      
    } else {
      text_add <- FALSE
      plotobj$group <- retrieved_by$value
    }
    
    if(is.null(legend_title)) legend_title <- by
    ggobj <- 
      ggplot(plotobj, aes(x=dim1, y=dim2, color=group)) +
      geom_point(alpha=list.ggplot$umap_alpha, size=list.ggplot$pointsize)
    
    if(class(retrieved_by$value)=="numeric"){
      
      ggobj <- 
        ggobj + 
        scale_color_viridis(name=legend_title)
      
    } else {
      
      ggobj <- 
        ggobj + 
        scale_color_manual(values=list.ggplot$colorblind_cols) +
        guides(colour=guide_legend(title=legend_title, override.aes=list(alpha=list.ggplot$umap_alpha, size=list.ggplot$legendsize)))
      
    }
    
    # add name of the by into the center of the by points:
    if(!is.null(label_size)){
      lsize <- label_size
    } else lsize <- list.ggplot$textsize
    
    if(text_add) {
      if(!text_use_label) text_fun <- geom_text else text_fun <- geom_label
      ggobj <- ggobj + text_fun(data=text_df, mapping=aes(x=x, y=y, label=group),
                                size=lsize, color=text_col)
    }
  }
  
  ggobj <- 
    ggobj +
    xlab(paste0(embed_onto, 1)) +
    ylab(paste0(embed_onto, 2))
  
  if(class(plotobj$values) != "numeric"){
    ggobj <- ggobj + theme(legend.position="top", legend.justification="center", 
                           legend.direction="horizontal")
  } 
  
  # either return plot+data or plot only:
  rm(sce)
  if(return_data){
    return(list(plot=ggobj, data=plotobj))
  } else {
    rm(plotobj)
    return(ggobj)
  }
  
}

#------------------------------------
# Custom function to plot scRNA-seq 
# QC metrics as boxplots with jittered
# points
#------------------------------------

plotQC <- function(data, x, y, x.lab=NULL, y.lab=NULL){
  suppressMessages(require(ggplot2))
  g <- ggplot(data, aes(x=get(x), y=get(y))) +
    geom_violin() + 
    geom_boxplot(width=0.1,outlier.size=0.05) +
    #geom_jitter(size=0.05, alpha=.1, height=0, color="grey", width=.25) +
    xlab(x.lab) + ylab(y.lab) +
    guides(x=guide_axis(angle=45))
  return(g)
}  

#------------------------------------
# Run gost() from g:profiler and 
# tidy upr results
#------------------------------------

run_gost <- function (InputGenes, BackgroundGenes=NULL, Species="mmusculus", 
                      Sources=c("KEGG", "REAC"), FDR.threshold=0.05, ...) {
  
  suppressMessages({
    require(gprofiler2)
    require(dplyr)
  })
  
  gg <- gost(query=InputGenes, custom_bg=BackgroundGenes, 
             organism=Species, exclude_iea=TRUE, evcodes=TRUE, 
             user_threshold=FDR.threshold, sources=Sources, ...)$result
  if (is.null(gg)) 
    return(NULL)
  gg <- data.frame(Term =gg$term_name, pvalue=gg$p_value, 
                   isize=gg$intersection_size, tsize=gg$term_size, Source=gg$source, 
                   Genes=gg$intersection)
  gg$Genes <- unlist(lapply(gg$Genes, 
                            function(x) paste(sort(strsplit(x, split=",")[[1]]), 
                                              collapse=",")))
  gg <- gg %>% arrange(pvalue) %>% filter(!grepl("KEGG|REACTOME", Term))
  return(gg)
}

#------------------------------------
# Select most variable rows, e.g. for PCA:
#------------------------------------

# select most variable rows based on rowVars:
#' @param counts a matrix of counts, should be log scale
#' @param ntop integer, how many genes to select
#' @param type either rv for simple rowVars() variance or mgv for modeling gene 
#' variance versus log2 counts with scran:
most_variable_rows <- function(counts, ntop=1000, type=c("rv", "mgv")){
  
  suppressMessages({
    require(matrixStats)
    require(scran)
  })
  
  invisible(match.arg(arg=class(counts)[1], choices=c("matrix")))
  invisible(match.arg(arg=class(ntop), choices=c("numeric", "integer")))
  type <- match.arg(type)
  
  if(type=="rv"){
    ## ntop most variable genes/regions:    
    rv <- matrixStats::rowVars(as.matrix(counts))
    
    if(length(rv) < ntop) message("Fewer than ntop genes in counts. Returning all genes.")
    
    selected <- head(order(rv, decreasing=TRUE), ntop)
    
    if(is.null(rownames(counts))){
      
      message("No rownames found in matrix, returning index of top variable rows")
      return(selected)
      
    } else return(rownames(counts)[selected])
  }
  
  if(type=="mgv"){
    
    if(is.null(rownames(counts))) {
      rwn=FALSE; message("No rownames found in matrix, returning index of top variable rows")
    }
    
    return(scran::getTopHVGs(scran::modelGeneVar(x=counts),n=ntop, row.names=rwn))
    
  }
  
}

#------------------------------------
# Opposite of intersect:
#------------------------------------

# https://www.r-bloggers.com/2011/11/outersect-the-opposite-of-rs-intersect-function/
outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}

#------------------------------------
# Wrapper around DESeq2::results() + DESeq2::lfcShrink to test and shrink a contrast,
# assumes DESeq() did already run on the dds:
#------------------------------------

#' results/lfcShrink wrapper

run_deseq2 <- function(dds, con, na_omit=TRUE, results_args=list(quiet=TRUE), lfcShrink_args=list(type="ashr", quiet=TRUE)){
  
  suppressMessages({
    require(ashr)
    require(DESeq2)
    require(tidyverse)
  })
  
  r <- do.call(DESeq2::results, c(list(object=dds, contrast=con), results_args))
  l <- do.call(DESeq2::lfcShrink(), c(list(dds=dds, contrast=con, res=r), lfcShrink_args)) %>%
    data.frame() %>%
    dplyr::mutate(baseMean=log2(baseMean+1))%>%
    data.frame(Gene=rownames(.), .)
  
  if(na_omit){
    return(na.omit(l))
  } else return(l)
  
}

#------------------------------------
# Customized version of fgsea::plotEnrichment
#------------------------------------

plot_fgsea <- function(pathway, stats, gseaParam=1, ticksSize=0.3,
                       lwd=0.5, color_line="darkorchid", color_segment="black",
                       xlab="ranked genes", ylab="enrichment score", add_cutoff=TRUE,
                       add_caption=FALSE){
  
  suppressMessages({
    require(fgsea)
    require(tidyverse)
  })
  
  l_pathway <- length(pathway)
  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
  statsAdj <- statsAdj/max(abs(statsAdj))
  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)
  gseaRes <- calcGseaStat(statsAdj, selectedStats=pathway, 
                          returnAllExtremes=TRUE)
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  diff <- (max(tops) - min(bottoms))/8
  x=y=NULL
  toPlot <- data.frame(x=c(0, xs, n + 1), y=c(0, ys, 0))
  
  minpos <- min(stats[stats>0])
  maxneg <- max(stats[stats<0])
  if(length(minpos) > 0) use_cutoff <- minpos else use_cutoff <- maxneg
  
  g <- 
    ggplot(toPlot, aes(x=x, y=y)) + 
    geom_hline(yintercept=0, color="black", lwd=lwd) +
    geom_line(color=color_line) + 
    geom_segment(data=data.frame(x=pathway), 
                 mapping=aes(x=x, y=-diff/2, xend=x, yend=diff/2), 
                 size=ticksSize, color=color_segment) + 
    labs(x=xlab, y=ylab)
  
  if(add_cutoff){
    g <- 
      g +
      geom_segment(data=data.frame(x=as.numeric(which(stats==use_cutoff)[1])), 
                   mapping=aes(x=x, y=-diff*1.1, xend=x, yend=diff*1.1), 
                   size=ticksSize*2)
  }
  g + 
    labs(caption=paste(l_pathway, "genes in pathway --", length(pathway), "found in ranking"))
  if(add_caption) {
    return(g + theme(plot.caption=element_text(hjust=0)))
  } else return(g)
}

#------------------------------------
# Take a vector and make all possible pairwise contrasts in edgeR-ish or DESeq2-ish style:
#------------------------------------

make_all_contrasts <- function (group, delim="_vs_", deseq2=FALSE, name="group"){
  
  suppressMessages(require(limma))
  
  group <- sort(unique(as.character(group)))
  if (sum(grepl("-", unlist(group))) > 0) 
    stop("There is a hyphen somewhere in group", call.=FALSE)
  cb <- combn(group, 2)
  Contrasts <- list()
  for (x in seq(1, ncol(cb))) Contrasts[[x]] <- paste0(cb[1, 
                                                          x], "-", cb[2, x])
  Contrasts <- limma::makeContrasts(contrasts=unlist(Contrasts), 
                                    levels=group)
  colnames(Contrasts) <- gsub("-", delim, colnames(Contrasts))
  message("Created ", ncol(Contrasts), " contrasts")
  if (deseq2) {
    return(sapply(colnames(Contrasts), function(x) {
      sp <- strsplit(x, split=delim)[[1]]
      return(c(name, sp[1], sp[2]))
    }, simplify=FALSE, USE.NAMES=TRUE))
  }
  else return(Contrasts)
}

#------------------------------------
# Convert ggplot sizes to gpar sizes
#------------------------------------

gg2gp <- function(x) ggplot2:::.pt*x

#------------------------------------
# Make identical elements unique but start at the first element, see
# https://stackoverflow.com/questions/7659891/r-make-unique-starting-in-1
#------------------------------------

make.unique.2 = function(x, sep='.'){
  ave(x, x, FUN=function(a){if(length(a) > 1){paste(a, 1:length(a), sep=sep)} else {a}})
}

#------------------------------------
# efficient matrixStats version of t(scale(t(x)))
#------------------------------------

rowScale <- function (x, center = TRUE, scale = TRUE, add_attr = FALSE, rows = NULL, 
                      cols = NULL, only_complete=TRUE){
  
  suppressMessages(require(matrixStats))
  
  if (!is.null(rows) && !is.null(cols)) {
    x <- x[rows, cols, drop = FALSE]
  }
  else if (!is.null(rows)) {
    x <- x[rows, , drop = FALSE]
  }
  else if (!is.null(cols)) {
    x <- x[, cols, drop = FALSE]
  }
  cm = rowMeans(x, na.rm = TRUE)
  if (scale) {
    csd = matrixStats::rowSds(x, center = cm)
  }
  else {
    csd = rep(1, length = length(cm))
  }
  if (!center) {
    cm = rep(0, length = length(cm))
  }
  x = (x - cm)/csd
  if (add_attr) {
    if (center) {
      attr(x, "scaled:center") <- cm
    }
    if (scale) {
      attr(x, "scaled:scale") <- csd
    }
  }
  
  if(only_complete) {
    return(x[complete.cases(x),])
  } else return(x)
}

#------------------------------------
# Split a string chr:1-10 at : and - and make a new df with chr-start-end with it:
#------------------------------------

string2chr <- function(x, pattern = ":|-",
                       chr.name = "chr",
                       start.name = "start",
                       end.name = "end"){
  
  require(data.table)
  
  splitted <- data.table::tstrsplit(x, pattern)
  if(length(splitted)!=3) stop("Split did not produce three elements")
  df<-data.frame(A=splitted[[1]], B=splitted[[2]], C=splitted[[3]])
  colnames(df) <- c(chr.name, start.name, end.name)
  df
  
}

#------------------------------------
# Scale a numeric matrix or df by quantiles
#------------------------------------

scale_by_quantile <- function (Counts, lower = 0, upper = 1){
  
  if (class(Counts)[1] != "matrix" & class(Counts)[1] != "data.frame") {
    stop("Counts must be a matrix or data.frame")
  }
  if (lower == 0 & upper == 1) 
    return(Counts)
  if (class(Counts)[1] == "data.frame") {
    cnames <- colnames(Counts)
    Counts <- as.matrix(Counts)
    colnames(Counts) <- cnames
  }
  if (upper < 1) {
    qt.upper <- as.numeric(quantile(Counts, upper, na.rm = TRUE))
    Counts[Counts > qt.upper] <- qt.upper
  }
  if (lower > 0) {
    qt.lower <- as.numeric(quantile(Counts, lower, na.rm = TRUE))
    Counts[Counts < qt.lower] <- qt.lower
  }
  return(Counts)
  
}

#------------------------------------
# Save plots as pdf and png
#------------------------------------

#' Save plots
#' 
#' Save plots stored in a list as pdf and png
#' 
#' @param plotlist a named list with ggplot objects or compatible grobs
#' @param outdir output directory name
#' @param prefix a prefix to give the plot names
#' @param overwrite logical, overwrite if plot already exists
#' @param width width of plot
#' @param height height of plot
#' @param res resolution for the png
#' @param units units for the png 
#' @param verbose logical whether to tell which plot it stored to which location
#' @param create_outdir logical, guess what
#' 
save_plots <- function(plotlist, outdir, prefix="",
                       overwrite=TRUE, width=7, height=7, res=300, units="in",
                       verbose=TRUE, create_outdir=TRUE, skip_png=FALSE){
  
  if(length(plotlist)==0){
    message("[Info] save_plots did not do anything because the input list is empty!")
    return(NULL)
  }
  
  if(is.null(names(plotlist))) 
    stop("plotlist has no names() -- names(plotlist) will be the filenames!")
  
  if(!dir.exists(outdir))
    if(create_outdir){
      message("Creating outdir")
      dir.create(outdir, recursive=TRUE)
    } else stop("Outdir does not exist and create_outdir=FALSE")
  
  for(p in names(plotlist)){
    
    tmp.plot    <- plotlist[[p]]
    tmp.name    <- paste0(outdir, "/", prefix, p)
    tmp.message <- gsub("//|///", "/", paste("Saving", p, "to", paste0(tmp.name, ".png/.pdf")))
    
    if(verbose) tmp.message
    
    # pdf:
    save_pdf <- paste0(tmp.name, ".pdf")
    if(!file.exists(save_pdf) | (file.exists(save_pdf) & overwrite)){
      pdf(paste0(tmp.name, ".pdf"), width=width, height=height)
      print(tmp.plot); dev.off()
    } else message(paste(save_pdf, "exists and overwrite is set to FALSE"))
    
    # png:
    if(!skip_png){
      save_png <- paste0(tmp.name, ".png")
      if(!file.exists(save_png) | (file.exists(save_png) & overwrite)){
        png(paste0(tmp.name, ".png"), width=width, height=height, res=res, units=units)
        print(tmp.plot); dev.off()
      } else message(paste(save_png, "exists and overwrite is set to FALSE"))
    }
    
    rm(tmp.plot, tmp.name, tmp.message)      
    
  }
  
}
