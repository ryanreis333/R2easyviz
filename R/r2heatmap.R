#' Generate a Heatmap of Marker Genes
#'
#' The `r2heatmap` function creates a heatmap of top marker genes from a Seurat object.
#' It visualizes the expression of selected markers across a downsampled set of cells,
#' with options for custom colors and grouping.
#'
#' @param seurat_obj A Seurat object.
#' @param FindAllMarkersObj A data frame from `FindAllMarkers` containing marker genes.
#' @param group_by A metadata column to group cells by in the heatmap. Default is `"celltype"`.
#' @param ncells An integer for the number of cells to sample per group. Default is 500.
#' @param viridis_color A logical indicating whether to use the Viridis color scale. Default is `TRUE`.
#' @param nfeatures An integer for the number of top marker genes per cluster to show. Default is 5.
#' @param arrange_by A column name in `FindAllMarkersObj` to order genes by.
#'   If `"dif"`, it's calculated as `pct.1 - pct.2`. Default is `"dif"`.
#' @param group_colors An optional named vector for custom group colors.
#'
#' @return A `ggplot` object representing the heatmap.
#' @export
#'
#' @examples
#' \dontrun{
#' # To create a dummy Seurat object for testing
#' pbmc_data <- matrix(rnorm(1000 * 20), nrow = 1000)
#' rownames(pbmc_data) <- paste0("Gene-", 1:1000)
#' colnames(pbmc_data) <- paste0("Cell-", 1:20)
#' seurat_obj <- CreateSeuratObject(counts = pbmc_data)
#' seurat_obj$celltype <- sample(c("T-cell", "B-cell"), 20, replace = TRUE)
#' Idents(seurat_obj) <- seurat_obj$celltype
#' markers_df <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#'
#' # Generate a heatmap
#' r2heatmap(seurat_obj, FindAllMarkersObj = markers_df, group_by = "celltype")
#' }
#'
#' @importFrom dplyr group_by slice_sample arrange desc mutate summarise pull
#' @importFrom rlang sym
#' @importFrom Seurat ScaleData DoHeatmap FetchData
#' @importFrom viridis viridis
#' @importFrom ggplot2 scale_fill_gradientn
#' @importFrom utils globalVariables
#'
r2heatmap <- function(seurat_obj,
                      FindAllMarkersObj,
                      group_by = "celltype",
                      ncells = 500,
                      viridis_color = TRUE,
                      nfeatures = 5,
                      arrange_by = "dif",
                      group_colors = NULL) {

  # Ensure FindAllMarkersObj is a data frame
  if (!is.data.frame(FindAllMarkersObj)) {
    stop("FindAllMarkersObj must be a data frame.")
  }

  # Extract metadata and add barcode column
  meta_data <- seurat_obj@meta.data
  meta_data$barcode <- rownames(meta_data)

  # Check if group_by exists
  if (!group_by %in% colnames(meta_data)) {
    stop(paste("Column", group_by, "not found in metadata."))
  }

  # Sample barcodes
  df_barcodes <- meta_data %>%
    group_by(!!sym(group_by)) %>%
    slice_sample(n = ncells, replace = TRUE)

  # Check and compute 'dif' if needed
  if (!(arrange_by %in% colnames(FindAllMarkersObj))) {
    if (arrange_by == "dif" && all(c("pct.1", "pct.2") %in% colnames(FindAllMarkersObj))) {
      FindAllMarkersObj <- FindAllMarkersObj %>%
        mutate(dif = pct.1 - pct.2)
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
  heatmap.markers <- unique(markers.subset$gene)

  # Scale data for the selected features
  seurat_obj <- ScaleData(seurat_obj, features = heatmap.markers)

  # Create heatmap
  plot <- DoHeatmap(
    object = seurat_obj,
    features = heatmap.markers,
    cells = df_barcodes$barcode,
    label = FALSE,
    group.by = group_by,
    group.colors = group_colors
  )

  # Apply optional Viridis color scale
  if (viridis_color) {
    plot <- plot + scale_fill_gradientn(colors = viridis(100))
  }

  return(plot)
}
