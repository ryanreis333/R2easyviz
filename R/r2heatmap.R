#' Generate a Heatmap of Marker Genes
#'
#' The `r2heatmap` function creates a heatmap of the top marker genes for a
#' given `Seurat` object. It allows you to visualize the expression of
#' selected marker genes across cells, with options for color scaling and
#' grouping. The function includes checks for column existence and
#' dynamically adjusts the number of cells sampled based on the available
#' data.
#'
#' @param seurat_obj A `Seurat` object containing the single-cell RNA-seq
#'   data. This object should include metadata with cell grouping
#'   information.
#' @param FindAllMarkersObj A data frame or tibble containing the results
#'   from `FindAllMarkers`, which includes information about marker genes.
#'   This should have columns for gene names, clusters, and the metric used
#'   for ranking.
#' @param group_by A character string specifying the metadata column used
#'   for grouping cells in the heatmap. Default is `"celltype"`.
#' @param ncells An integer specifying the number of cells to sample per
#'   group for the heatmap. If `ncells` is greater than the minimum number
#'   of cells in any group, it will be adjusted to the minimum. Default is
#'   `500`.
#' @param viridis_color A logical value indicating whether to use the
#'   Viridis color scale for the heatmap. If `TRUE`, the Viridis color scale
#'   will be applied; otherwise, a default color scale will be used. Default
#'   is `TRUE`.
#' @param nfeatures An integer specifying the number of top marker genes to
#'   include in the heatmap per cluster. Default is `5`.
#' @param arrange_by A character string specifying the column name in
#'   `FindAllMarkersObj` used for ordering marker genes. If `"dif"`, and the
#'   `dif` column is not present, it will be calculated as the difference
#'   between `pct.1` and `pct.2`. Default is `"dif"`.
#' @param barcode_column A character string specifying the column name in
#'   the metadata that contains barcode information. Default is `"barcodes"`.
#' @param group_colors An optional named vector specifying manual colors for
#'   `group_by` groups. Names must match the group names exactly. Default is
#'   `NULL`.
#'
#' @return A `ggplot` object representing the heatmap of the selected marker
#'   genes.
#' @export
#'
#' @examples
#' \dontrun{
#' # Top 5 marker genes, Viridis color scale
#' heatmap_plot <- r2heatmap(
#'   seurat_obj = pbmc,
#'   FindAllMarkersObj = markers_df,
#'   group_by = "celltype",
#'   nfeatures = 5
#' )
#'
#' # Top 10 marker genes, default color scale
#' heatmap_plot <- r2heatmap(
#'   seurat_obj = pbmc,
#'   FindAllMarkersObj = markers_df,
#'   group_by = "celltype",
#'   nfeatures = 10,
#'   viridis_color = FALSE
#' )
#' }
#'
#' @importFrom dplyr group_by slice_sample arrange desc mutate slice_head
#' @importFrom rlang sym
#' @importFrom Seurat ScaleData DoHeatmap
#' @importFrom ggplot2 scale_fill_gradientn
#' @importFrom viridis viridis
r2heatmap <- function(seurat_obj,
                      FindAllMarkersObj,
                      group_by = "celltype",
                      ncells = 500,
                      viridis_color = TRUE,
                      nfeatures = 5,
                      arrange_by = "dif",
                      barcode_column = "barcodes",
                      group_colors = NULL) {

  # Ensure FindAllMarkersObj is a data frame
  if (!is.data.frame(FindAllMarkersObj)) {
    stop("FindAllMarkersObj must be a data frame.")
  }

  # Extract metadata and sample barcodes
  meta_data <- seurat_obj[[]]
  df_barcodes <- meta_data %>%
    dplyr::group_by(!!rlang::sym(group_by)) %>%
    dplyr::slice_sample(n = ncells, replace = FALSE)

  # Check if barcode_column exists
  if (!(barcode_column %in% colnames(df_barcodes))) {
    warning(paste("Column", barcode_column, "not found in metadata."))
    return(NULL)
  }

  # Check if arrange_by exists and create 'dif' if needed
  if (!(arrange_by %in% colnames(FindAllMarkersObj))) {
    if (arrange_by == "dif" && all(c("pct.1", "pct.2") %in% colnames(FindAllMarkersObj))) {
      FindAllMarkersObj <- FindAllMarkersObj %>%
        dplyr::mutate(dif = .data$pct.1 - .data$pct.2)
      arrange_by <- "dif"
    } else {
      stop(paste("Column", arrange_by, "not found in 'FindAllMarkersObj' and cannot calculate 'dif' column."))
    }
  }

  # Subset marker genes
  markers.subset <- FindAllMarkersObj %>%
    dplyr::group_by(.data$cluster) %>%
    dplyr::arrange(dplyr::desc(!!rlang::sym(arrange_by))) %>%
    dplyr::slice_head(n = nfeatures)

  # Extract features for the heatmap
  heatmap.markers <- markers.subset$gene

  # Scale data for the selected features
  seurat_obj <- Seurat::ScaleData(seurat_obj, features = heatmap.markers)

  # Create heatmap
  plot <- Seurat::DoHeatmap(
    object = seurat_obj,
    features = heatmap.markers,
    cells = df_barcodes[[barcode_column]],
    label = FALSE,
    group.by = group_by,
    group.colors = group_colors
  )

  # Apply optional Viridis color scale for expression values
  if (viridis_color) {
    plot <- plot + ggplot2::scale_fill_gradientn(colors = viridis::viridis(100))
  }

  return(plot)
}
