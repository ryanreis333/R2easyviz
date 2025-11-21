#' Calculate Cell Type Proportions from Seurat Object
#'
#' This function accepts a Seurat object, extracts the metadata, and calculates the
#' proportion of each cell type/cluster within each sample. It explicitly handles
#' missing cell types by filling them with 0.
#'
#' @param seurat_obj A Seurat object containing single-cell data.
#' @param sample_col String. The column name in \code{meta.data} representing the sample (e.g., "orig.ident").
#' @param celltype_col String. The column name in \code{meta.data} representing the cell type or cluster (e.g., "seurat_clusters").
#' @param group_col String (Optional). A column in \code{meta.data} to keep for grouping (e.g., "condition").
#'
#' @return A dataframe in long format with columns:
#' \itemize{
#'   \item The sample column
#'   \item The cell type column
#'   \item \code{n}: Raw count
#'   \item \code{proportion}: Relative frequency (0-1)
#'   \item The grouping column (if provided)
#' }
#' @import dplyr
#' @import tidyr
#' @import Seurat
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming 'pbmc' is your Seurat object
#' props <- r2prop_df(
#'   seurat_obj = pbmc,
#'   sample_col = "orig.ident",
#'   celltype_col = "seurat_clusters",
#'   group_col = "treatment"
#' )
#' }
r2prop_df <- function(seurat_obj, sample_col, celltype_col, group_col = NULL) {

  # 1. Input Validation
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Error: The input 'seurat_obj' must be a Seurat object.")
  }

  # Check if required packages are installed
  if (!requireNamespace("dplyr", quietly = TRUE) || !requireNamespace("tidyr", quietly = TRUE)) {
    stop("The packages 'dplyr' and 'tidyr' are required for this function.")
  }

  # 2. Extract Metadata
  # We use [[]] or @meta.data to get the dataframe
  meta_df <- seurat_obj@meta.data

  # Verify columns exist
  if (!sample_col %in% colnames(meta_df)) {
    stop(paste("Column", sample_col, "not found in Seurat metadata."))
  }
  if (!celltype_col %in% colnames(meta_df)) {
    stop(paste("Column", celltype_col, "not found in Seurat metadata."))
  }

  # 3. Calculate Counts
  # Group by sample and celltype to get raw counts
  counts_df <- meta_df %>%
    dplyr::group_by(.data[[sample_col]], .data[[celltype_col]]) %>%
    dplyr::summarise(n = dplyr::n(), .groups = 'drop')

  # 4. Fill Zeros (Crucial Step)
  # Use tidyr::complete to ensure samples with 0 cells of a specific type are recorded as 0
  full_df <- counts_df %>%
    tidyr::complete(.data[[sample_col]], .data[[celltype_col]], fill = list(n = 0))

  # 5. Calculate Proportions
  prop_df <- full_df %>%
    dplyr::group_by(.data[[sample_col]]) %>%
    dplyr::mutate(
      total_cells = sum(n),
      proportion = n / total_cells
    ) %>%
    dplyr::ungroup()

  # 6. Re-attach Group Metadata (Optional)
  if (!is.null(group_col)) {
    if (!group_col %in% colnames(meta_df)) {
      stop(paste("Column", group_col, "not found in Seurat metadata."))
    }

    # Create a unique mapping of Sample -> Group
    # We assume a sample ID corresponds to a single condition
    meta_map <- meta_df %>%
      dplyr::select(dplyr::all_of(c(sample_col, group_col))) %>%
      dplyr::distinct()

    # Check for mapping issues (e.g., one sample ID having multiple conditions)
    if (any(duplicated(meta_map[[sample_col]]))) {
      warning(paste("The sample column", sample_col, "is not unique with respect to", group_col,
                    ". Some samples may have mixed conditions. Duplicate rows will be created."))
    }

    prop_df <- prop_df %>%
      dplyr::left_join(meta_map, by = sample_col)
  }

  return(prop_df)
}
