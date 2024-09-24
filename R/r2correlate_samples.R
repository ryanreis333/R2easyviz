#' Correlate Aggregated Expression Data from a Seurat Object
#'
#' This function aggregates gene expression data from a Seurat object based on a specified grouping factor
#' (e.g., sample or cluster) and computes the correlation matrix between the groups. Optionally, it can return
#' a heatmap visualization of the correlation matrix.
#'
#' @param seurat_obj A Seurat object containing the single-cell RNA-seq data.
#' @param group_by A character string specifying the metadata column in the Seurat object to group cells by.
#'   This can be a column like `"orig.ident"` or any other metadata column in `seurat_obj@meta.data`.
#' @param assay A character string indicating the assay to use from the Seurat object. Default is `"RNA"`.
#' @param return_heatmap A logical value. If `TRUE`, the function returns a heatmap of the correlation matrix using the `pheatmap` package.
#'   If `FALSE`, the function returns the correlation matrix itself. Default is `TRUE`.
#'
#' @return If `return_heatmap = TRUE`, returns a heatmap of the correlation matrix. If `return_heatmap = TRUE`,
#'   returns a matrix where each element represents the correlation between two groups' aggregated gene expression.
#'
#' @details The function aggregates gene expression data based on the `group_by` parameter using the `AggregateExpression()`
#' function from Seurat. It then computes the correlation matrix between groups using Pearson's correlation
#' (`cor` function with `use = "complete.obs"` to handle missing values). If requested, the function will visualize the
#' correlation matrix as a heatmap using the `pheatmap` package.
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' pbmc <- CreateSeuratObject(pbmc_data)
#' correlation_matrix <- r2correlate_samples(seurat_obj = pbmc, group_by = "orig.ident", assay = "RNA", return_heatmap = FALSE)
#'
#' # Return a heatmap instead of the matrix
#' r2correlate_samples(seurat_obj = pbmc, group_by = "orig.ident", return_heatmap = TRUE)
#' }
#'
#' @import Seurat pheatmap
#' @export
r2correlate_samples <- function(seurat_obj, group_by = "orig.ident", assay = "RNA", return_heatmap = T){

  # Aggregate expression based on group_by and assay
  agg.seurat <- AggregateExpression(object = seurat_obj,
                                    group.by = group_by,
                                    return.seurat = TRUE,
                                    assays = assay)

  # Assign column names based on the group_by metadata
  colnames(agg.seurat) <- unique(seurat_obj[[group_by]][,1])

  # Extract the aggregated expression matrix
  mat <- GetAssayData(agg.seurat, slot = "data")

  # Initialize a correlation matrix
  cor_matrix <- matrix(nrow = ncol(mat), ncol = ncol(mat))

  # Compute correlation matrix
  for(i in 1:ncol(mat)){
    for (j in 1:ncol(mat)) {
      col1 <- mat[, i]
      col2 <- mat[, j]
      cor_matrix[i, j] <- cor(col1, col2, use = "complete.obs")
    }
  }

  # Add row and column names to the correlation matrix
  colnames(cor_matrix) <- colnames(agg.seurat)
  rownames(cor_matrix) <- colnames(agg.seurat)

  # Return heatmap if requested
  if(return_heatmap){
    return(pheatmap::pheatmap(cor_matrix, show_rownames = TRUE, show_colnames = TRUE))
  }

  # Otherwise return the correlation matrix
  return(cor_matrix)
}
