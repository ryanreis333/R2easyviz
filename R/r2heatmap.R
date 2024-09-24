#' Generate a Heatmap of Marker Genes
#'
#' The `r2heatmap` function creates a heatmap of the top marker genes for a given `Seurat` object. It allows you to visualize the expression of selected marker genes across cells, with options for color scaling and grouping. The function includes checks for column existence and dynamically adjusts the number of cells sampled based on the available data.
#'
#' @param seurat_obj A `Seurat` object containing the single-cell RNA-seq data. This object should include metadata with cell grouping information.
#' @param FindAllMarkersObj A data frame or tibble containing the results from `FindAllMarkers`, which includes information about marker genes. This should have columns for gene names, clusters, and the metric used for ranking.
#' @param group_by A character string specifying the metadata column used for grouping cells in the heatmap. Default is `"celltype"`.
#' @param ncells An integer specifying the number of cells to sample per group for the heatmap. If `ncells` is greater than the minimum number of cells in any group, it will be adjusted to the minimum. Default is `500`.
#' @param viridis_color A logical value indicating whether to use the Viridis color scale for the heatmap. If `TRUE`, the Viridis color scale will be applied; otherwise, a default color scale will be used. Default is `TRUE`.
#' @param nfeatures An integer specifying the number of top marker genes to include in the heatmap per cluster. Default is `5`.
#' @param arrange_by A character string specifying the column name in `FindAllMarkersObj` used for ordering marker genes. If `"dif"`, and the `dif` column is not present, it will be calculated as the difference between `pct.1` and `pct.2`. Default is `"dif"`.
#' @param barcode_column A character string specifying the column name in the metadata that contains barcode information. Default is `"barcodes"`.
#'
#' @return A `ggplot` object representing the heatmap of the selected marker genes.
#' @export
#'
#' @examples
#' # Generate a heatmap for the top 5 marker genes, using the Viridis color scale
#' heatmap_plot <- r2heatmap(seurat_obj = pbmc, FindAllMarkersObj = markers_df, group_by = "celltype", nfeatures = 5)
#'
#' # Generate a heatmap for the top 10 marker genes, without using the Viridis color scale
#' heatmap_plot <- r2heatmap(seurat_obj = pbmc, FindAllMarkersObj = markers_df, group_by = "celltype", nfeatures = 10, viridis_color = FALSE)
#'
#' @importFrom dplyr group_by slice_sample arrange desc mutate summarise pull
#' @importFrom rlang sym
#' @importFrom Seurat ScaleData DoHeatmap
#' @importFrom viridis viridis
#'
r2heatmap <- function(seurat_obj, FindAllMarkersObj, group_by = "celltype", ncells = 500, viridis_color = TRUE, nfeatures = 5, arrange_by = "dif", barcode_column = "barcodes") {

  # Ensure FindAllMarkersObj is a data frame
  if (!is.data.frame(FindAllMarkersObj)) {
    stop("FindAllMarkersObj must be a data frame.")
  }

  # Extract metadata and sample barcodes
  meta_data <- seurat_obj[[]]
  df_barcodes <- meta_data %>%
    group_by(!!sym(group_by)) %>%
    slice_sample(n = ncells, replace = FALSE)

  # Check if barcode_column exists
  if (!(barcode_column %in% colnames(df_barcodes))) {
    warning(paste("Column", barcode_column, "not found in metadata."))
    return(NULL)  # Exit the function if barcode_column is not found
  }

  # Check if arrange_by exists and create 'dif' if needed
  if (!(arrange_by %in% colnames(FindAllMarkersObj))) {
    if (arrange_by == "dif" && all(c("pct.1", "pct.2") %in% colnames(FindAllMarkersObj))) {
      FindAllMarkersObj <- FindAllMarkersObj %>%
        mutate(dif = pct.1 - pct.2)
      arrange_by <- "dif"  # Ensure we're using the newly created 'dif' column
    } else {
      stop(paste("Column", arrange_by, "not found in 'FindAllMarkersObj' and cannot calculate 'dif' column."))
    }
  }

  # Subset marker genes
  markers.subset <- FindAllMarkersObj %>%
    group_by(cluster) %>%
    arrange(desc(!!sym(arrange_by))) %>%
    slice_head(n = nfeatures)

  # Extract features for the heatmap
  heatmap.markers <- markers.subset$gene

  # Scale data for the selected features
  seurat_obj <- ScaleData(seurat_obj, features = heatmap.markers)

  # Create heatmap with optional Viridis color scale
  if (viridis_color) {
    plot <- DoHeatmap(object = seurat_obj, features = heatmap.markers, cells = df_barcodes[[barcode_column]], label = FALSE, group.by = group_by) +
      scale_fill_gradientn(colors = viridis(100))
  } else {
    plot <- DoHeatmap(object = seurat_obj, features = heatmap.markers, cells = df_barcodes[[barcode_column]], label = FALSE, group.by = group_by)
  }

  return(plot)
}
