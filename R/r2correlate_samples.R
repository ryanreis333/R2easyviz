#' Correlate Aggregated Expression Data from a Seurat Object
#'
#' This function aggregates gene expression data from a Seurat object based on a specified grouping factor
#' (e.g., sample or cluster) and computes a correlation matrix between the groups.
#' The function can return either the correlation matrix or a heatmap visualization.
#'
#' @param seurat_obj A Seurat object.
#' @param group_by A character string specifying the metadata column to group cells by.
#'   Must be a valid column name in the Seurat object's metadata.
#' @param assay A character string indicating which assay to use. Default is `"RNA"`.
#' @param return_heatmap A logical value. If `TRUE` (default), a heatmap of the
#'   correlation matrix is returned. If `FALSE`, the correlation matrix is returned.
#'
#' @return A `pheatmap` object if `return_heatmap` is `TRUE`, or a correlation matrix if `FALSE`.
#'
#' @details
#' This function first aggregates gene expression data for each group specified by the `group_by`
#' parameter using `Seurat::AggregateExpression`. Then, it calculates the Pearson correlation
#' between the aggregated expression profiles of the groups.
#'
#' @examples
#' \dontrun{
#' # To create a dummy Seurat object for testing
#' pbmc_data <- matrix(rnorm(1000 * 20), nrow = 1000)
#' rownames(pbmc_data) <- paste0("Gene-", 1:1000)
#' colnames(pbmc_data) <- paste0("Cell-", 1:20)
#' seurat_obj <- CreateSeuratObject(counts = pbmc_data)
#' seurat_obj$orig.ident <- sample(c("SampleA", "SampleB"), 20, replace = TRUE)
#'
#' # Get correlation matrix
#' cor_matrix <- r2correlate_samples(seurat_obj, group_by = "orig.ident", return_heatmap = FALSE)
#'
#' # Get heatmap
#' r2correlate_samples(seurat_obj, group_by = "orig.ident")
#' }
#'
#' @importFrom Seurat AggregateExpression GetAssayData FetchData
#' @importFrom pheatmap pheatmap
#' @importFrom stats cor
#' @export
r2correlate_samples <- function(seurat_obj, group_by = "orig.ident", assay = "RNA", return_heatmap = TRUE) {

  # Check if group_by exists in the metadata
  if (!group_by %in% colnames(seurat_obj[[]])) {
    stop("group_by column '", group_by, "' not found in Seurat object metadata.")
  }

  # Aggregate expression
  agg_seurat <- AggregateExpression(
    object = seurat_obj,
    group.by = group_by,
    assays = assay,
    return.seurat = TRUE
  )

  # Get aggregated expression matrix
  mat <- GetAssayData(agg_seurat, assay = assay, layer = "data")

  # Calculate correlation matrix
  cor_matrix <- cor(mat, use = "complete.obs")

  # Return heatmap or matrix
  if (return_heatmap) {
    return(pheatmap::pheatmap(cor_matrix, show_rownames = TRUE, show_colnames = TRUE))
  } else {
    return(cor_matrix)
  }
}
