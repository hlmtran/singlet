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
                               max.terms.per.factor = 3,
                               dropcommon = TRUE,
                               gsea.name = "gsea",
                               truncate.terms = 48, ...) {
  df <- object@reductions[[reduction]]@misc[[gsea.name]][["padj"]]
  
  p <- GSEAHeatmap.matrix(
    df,
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
                                             max.terms.per.factor = 3,
                                             dropcommon = TRUE,
                                             gsea.name = "gsea",
                                             truncate.terms = 48, ...) {
  df <- object@metadata[[reduction]][[gsea.name]][["padj"]]
  
  p <- GSEAHeatmap.matrix(
    df,
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
                            max.terms.per.factor = 3,
                            dropcommon = TRUE,
                            gsea.name = "gsea",
                            truncate.terms = 48, ...) {
  df <- object@misc[[gsea.name]][["padj"]]
  
  p <- GSEAHeatmap.matrix(
    df,
    max.terms.per.factor = max.terms.per.factor,
    dropcommon = dropcommon,
    truncate.terms = truncate.terms
  )
  
  return(p)
}

.S3method("GSEAHeatmap", "nmf", GSEAHeatmap.nmf)


#' @rdname GSEAHeatmap
#' @method GSEAHeatmap matrix
#' @export
GSEAHeatmap.matrix <- function(object,
                               max.terms.per.factor = 3,
                               dropcommon = TRUE,
                               truncate.terms = 48, ...) {
  df <- selectGSEATerms(
    object,
    max.terms.per.factor = max.terms.per.factor,
    dropcommon = dropcommon,
    truncate.terms = truncate.terms
  )

  df <- reshape2::melt(df)
  p <- plotGSEAHeatmap(df)
  
  return(p)
}

.S3method("GSEAHeatmap", "matrix", GSEAHeatmap.matrix)



# GSEAHeatmap <- function(object, reduction = "nmf", max.terms.per.factor = 3, dropcommon = TRUE,gsea.name = "gsea",truncate.terms=48) {
# #TODO: need to check if gsea.name is in object
#   if (is(object, "Seurat")) {
#     df <- object@reductions[[reduction]]@misc[[gsea.name]][["padj"]]
#   } else if (is(object, "nmf")) {
#     df <- object@misc[[gsea.name]][["padj"]]
#   }
#   
#   
#   df = selectGSEATerms(df,max.terms.per.factor,dropcommon,truncate.terms)
#   df <- reshape2::melt(df)
#   p <- plotGSEAHeatmap(df)
# 
#   return(p) 
# 
# }

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

selectGSEATerms = function(df, max.terms.per.factor = 3,dropcommon = TRUE,truncate.terms=48){
  df = df[(rowSums(df,na.rm=T) > 0 ),]
  # markers for each factor based on the proportion of signal in that factor
  df2 <- as.matrix(Diagonal(x = 1 / rowSums(df,na.rm = TRUE)) %*% df)
  
  # see https://github.com/zdebruine/singlet/issues/26
  # thanks to @earbebarnes
  rownames(df2) <- rownames(df) #add row names to df2
  
  terms <- c()
  for (i in 1:ncol(df2)) {
    terms_i <- df[, i]
    idx <- terms_i > -log10(0.05)
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
  df <- df[terms, ]
  #could be earlier, but this lets us exclude NAs where 0 could be < max.terms 
  df[is.na(df)] = 0
  
  if(!is.null(truncate.terms)){
    rownames(df) <- sapply(rownames(df), function(x) {
      ifelse(nchar(x) > truncate.terms, paste0(substr(x, 1, truncate.terms), "..."), x)
    })
  }
  
  if (dropcommon) { 
    # remove terms that are broadly significant
    v <- which((rowSums(df > -log10(0.05)) > (ncol(df) / 2)))
    if (length(v) > 0) df <- df[-v, ]
  }
  
  return(df)
}