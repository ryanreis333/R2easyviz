#' Downsample Seurat Object by Group
#'
#' `r2sample_seurat` is designed to downsample a Seurat object by randomly sampling
#' a specified number of cells (`n`) from each group defined by a metadata variable (`group_var`).
#' The function returns a new Seurat object containing only the sampled cells.
#'
#' @param seurat_obj A Seurat object containing the single-cell RNA-seq data.
#' @param group_var A character string representing the metadata column used to define the groups for sampling. Default is `"orig.ident"`.
#' @param n An integer specifying the number of cells to sample from each group. Default is 500.
#' @param barcode_column A character string representing the column in the Seurat metadata that contains the cell barcodes. Default is `"barcodes"`.
#'
#' @return A downsampled Seurat object containing the specified number of cells from each group.
#' @export
#'
#' @examples
#' # Downsample 100 cells from each sample group defined by `orig.ident`
#' pbmc.sample <- r2_sample_seurat(seurat_obj = pbmc, group_var = "orig.ident", n = 100)
#'
#' # Downsample 200 cells using a custom grouping variable and barcode column
#' seurat_downsampled <- r2sample_seurat(seurat_obj = my_seurat, group_var = "sample_id", n = 200, barcode_column = "cell_ids")
#'
r2sample_seurat <- function(seurat_obj, group_var = "orig.ident", n = 500, barcode_column = "barcodes") {
  # Extract the metadata from the Seurat object
  meta_data <- seurat_obj[[]]

  # Check if the specified barcode column exists in the metadata
  if (!(barcode_column %in% colnames(meta_data))) {
    stop(paste("Barcode column", barcode_column, "not found in Seurat object metadata."))
  }

  # Sample 'n' cells from each group specified in the group_var
  df_barcodes <- meta_data %>%
    group_by(!!sym(group_var)) %>%
    slice_sample(n = n)

  # Display the count of cells sampled from each group
  print(table(df_barcodes[[group_var]]))

  # Subset the Seurat object based on the sampled barcodes
  seurat_subset <- subset(seurat_obj, cells = df_barcodes[[barcode_column]])

  # Return the downsampled Seurat object
  return(seurat_subset)
}
