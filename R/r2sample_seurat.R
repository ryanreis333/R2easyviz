#' Downsample Seurat Object by Group
#'
#' `r2sample_seurat` is designed to downsample a Seurat object by randomly sampling
#' a specified number of cells (`n`) from each group defined by a metadata variable (`group_var`).
#' The function returns a new Seurat object containing only the sampled cells.
#'
#' @param seurat_obj A Seurat object containing the single-cell RNA-seq data.
#' @param group_var A character string representing the metadata column used to define the groups for sampling. Default is `"orig.ident"`.
#' @param n An integer specifying the number of cells to sample from each group. Default is 500.
#'
#' @return A downsampled Seurat object containing the specified number of cells from each group.
#' @export
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
#' # Downsample 10 cells from each sample group defined by `orig.ident`
#' pbmc_sampled <- r2sample_seurat(seurat_obj = seurat_obj, group_var = "orig.ident", n = 10)
#' }
#'
#' @importFrom Seurat FetchData
#' @importFrom dplyr %>% group_by slice_sample
#' @importFrom rlang sym
#'
r2sample_seurat <- function(seurat_obj, group_var = "orig.ident", n = 500) {
  # Extract the metadata from the Seurat object
  meta_data <- seurat_obj@meta.data
  meta_data$barcode <- rownames(meta_data)


  # Check if the specified group_var column exists in the metadata
  if (!(group_var %in% colnames(meta_data))) {
    stop(paste("Grouping variable", group_var, "not found in Seurat object metadata."))
  }

  # Sample 'n' cells from each group specified in the group_var
  df_barcodes <- meta_data %>%
    group_by(!!sym(group_var)) %>%
    slice_sample(n = n, replace = TRUE) # use replace = TRUE to avoid errors if n > number of cells in a group

  # Display the count of cells sampled from each group
  print(table(df_barcodes[[group_var]]))

  # Subset the Seurat object based on the sampled barcodes
  seurat_subset <- subset(seurat_obj, cells = df_barcodes$barcode)

  # Return the downsampled Seurat object
  return(seurat_subset)
}
