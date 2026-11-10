#' Plot GSEA results on a heatmap
#'
#' Plot top GSEA terms for each NMF factor on a heatmap
#'
#' @param object Seurat or RcppML::nmf object
#' @param reduction a dimensional reduction for which GSEA analysis has been performed
#' @param max.terms.per.factor show this number of top terms for each factor
#' @param dropcommon  drop broadly enriched terms across factors? (TRUE) 
#' @param gsea.name gsea name, gsea by default
#' @param truncate.terms max number of characters to keep in a gene set term. NULL to keep all.
#'
#' @return ggplot2 object
#'
#' @export
#'
GSEAHeatmap <- function(object, ...) {
  UseMethod("GSEAHeatmap")
}

#' @rdname GSEAHeatmap
#' @method GSEAHeatmap Seurat
#' @export
GSEAHeatmap.Seurat <- function(object, reduction = "nmf",
                               collapse = F,
                               padj.sig = 0.05,
                               max.terms.per.factor = 3,
                               dropcommon = TRUE,
                               gsea.name = "gsea",
                               truncate.terms = 48, ...) {
  gsea <- object@reductions[[reduction]]@misc[[gsea.name]]
  
  
  
  p <- GSEAHeatmap(
    object = gsea,
    collapse = collapse,
    padj.sig = padj.sig,
    max.terms.per.factor = max.terms.per.factor,
    dropcommon = dropcommon,
    truncate.terms = truncate.terms
  )
  
  return(p)
}
.S3method("GSEAHeatmap", "Seurat", GSEAHeatmap.Seurat)

#' @rdname GSEAHeatmap
#' @method GSEAHeatmap SingleCellExperiment
#' @export
GSEAHeatmap.SingleCellExperiment <- function(object, reduction = "nmf",
                                             collapse = F,
                                             padj.sig = 0.05,
                                             max.terms.per.factor = 3,
                                             dropcommon = TRUE,
                                             gsea.name = "gsea",
                                             truncate.terms = 48, ...) {
  gsea <- object@metadata[[reduction]][[gsea.name]]
  
  p <- GSEAHeatmap(
    object = gsea,
    collapse = collapse,
    padj.sig = padj.sig,
    max.terms.per.factor = max.terms.per.factor,
    dropcommon = dropcommon,
    truncate.terms = truncate.terms
  )
  
  return(p)
}

.S3method("GSEAHeatmap", "SingleCellExperiment", GSEAHeatmap.SingleCellExperiment)

#' @rdname GSEAHeatmap
#' @method GSEAHeatmap nmf
#' @export
GSEAHeatmap.nmf <- function(object,
                            collapse = F,
                            padj.sig = 0.05,
                            max.terms.per.factor = 3,
                            dropcommon = TRUE,
                            gsea.name = "gsea",
                            truncate.terms = 48, ...) {
  gsea <- object@misc[[gsea.name]]
  
  p <- GSEAHeatmap(
    object = gsea,
    collapse = collapse,
    padj.sig = padj.sig,
    max.terms.per.factor = max.terms.per.factor,
    dropcommon = dropcommon,
    truncate.terms = truncate.terms
  )
  
  return(p)
}

.S3method("GSEAHeatmap", "nmf", GSEAHeatmap.nmf)

#' @rdname GSEAHeatmap
#' @method GSEAHeatmap list
#' @importFrom fgsea collapsePathways
#' @export
GSEAHeatmap.list <- function(object,
                             collapse = F,
                             padj.sig = 0.05,
                               max.terms.per.factor = 3,
                               dropcommon = TRUE,
                               truncate.terms = 48, ...) {
  
  if(collapse){
    collapsed =  collapseGSEA(object, padj.sig = padj.sig)
    main_pathways = getMainPathways(collapsed)
    formatted =  formatGSEA(
      object,
      pathways = main_pathways
    )
  } else{
    formatted = formatGSEA(object)
  }
  
  p <- GSEAHeatmap.matrix(
    object = formatted[["padj"]],
    padj.sig = padj.sig,
    max.terms.per.factor = max.terms.per.factor,
    dropcommon = dropcommon,
    truncate.terms = truncate.terms
  )
  
  return(p)
}

.S3method("GSEAHeatmap", "list", GSEAHeatmap.list)



#' @rdname GSEAHeatmap
#' @method GSEAHeatmap matrix
#' @export
GSEAHeatmap.matrix <- function(object,
                               padj.sig = 0.05,
                               max.terms.per.factor = 3,
                               dropcommon = TRUE,
                               truncate.terms = 48, ...) {
  df <- selectGSEATerms(
    df = object,
    padj.sig=padj.sig,
    max.terms.per.factor = max.terms.per.factor,
    dropcommon = dropcommon,
    truncate.terms = truncate.terms
  )
  
  if (nrow(df) == 0) {
    stop("No significant GSEA pathways remain after filtering")
  }
  
  df <- reshape2::melt(df)
  p <- plotGSEAHeatmap(df)
  
  return(p)
}

.S3method("GSEAHeatmap", "matrix", GSEAHeatmap.matrix)


plotGSEAHeatmap = function(df){
  p <- ggplot(df, aes(Var2, Var1, fill = value)) +
    geom_tile() +
    scale_fill_viridis_c(option = "B") +
    theme_classic() +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(
      x = "NMF factor",
      y = "Gene Set",
      fill = "FDR\n(-log10)"
    ) +
    theme(
      axis.text.y = element_text(size = 6),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    ) + 
    NULL 
  
  return(p) 
}

selectGSEATerms = function(df, padj.sig = 0.05 ,max.terms.per.factor = 3,dropcommon = TRUE,truncate.terms=48){
  df = df[(rowSums(df,na.rm=T) > 0 ), , drop = FALSE]
  # markers for each factor based on the proportion of signal in that factor
  df2 <- as.matrix(Diagonal(x = 1 / rowSums(df,na.rm = TRUE)) %*% df)
  
  # see https://github.com/zdebruine/singlet/issues/26
  # thanks to @earbebarnes
  rownames(df2) <- rownames(df) #add row names to df2
  
  terms <- c()
  for (i in 1:ncol(df2)) {
    terms_i <- df[, i]
    idx <- terms_i > -log10(padj.sig)
    terms_i <- terms_i[idx]
    terms_j <- df2[idx, i]
    v <- sort(terms_j, decreasing = TRUE)
    v <- v[!is.na(v)] #making this explicit
    if (length(v) > max.terms.per.factor) {
      terms <- c(terms, names(v)[1:max.terms.per.factor])
    } else {
      terms <- c(terms, names(v))
    }
  }
  terms <- unique(terms)
  df <- df[terms, , drop = FALSE]
  #could be earlier, but this lets us exclude NAs where 0 could be < max.terms 
  df[is.na(df)] = 0
  
  if(!is.null(truncate.terms)){
    rownames(df) <- sapply(rownames(df), function(x) {
      ifelse(nchar(x) > truncate.terms, paste0(substr(x, 1, truncate.terms), "..."), x)
    })
  }
  
  if (dropcommon) { 
    # remove terms that are broadly significant
    v <- which((rowSums(df > -log10(padj.sig)) > (ncol(df) / 2)))
    if (length(v) > 0) df <- df[-v, , drop = FALSE]
  }
  
  return(df)
}