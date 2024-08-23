#' Proportional Plot of Cell Types Across Samples
#'
#' This function generates a bar plot showing the proportion of each cell type or specified metadata category across different samples or conditions from a Seurat object. By default, it uses the identities of the cells and the "orig.ident" column for grouping and splitting.
#'
#' @param seurat_object A Seurat object containing the data to be plotted. The object should include metadata columns corresponding to the cell types and sample identifiers.
#' @param celltype A string specifying the column name in the Seurat object's metadata that contains the cell type or other categorical metadata to be plotted. By default, it uses `Idents(seurat_object)` which refers to the identities of the cells.
#' @param orig.ident A string specifying the column name in the Seurat object's metadata that contains the sample identifiers or conditions by which the data should be split. The default is "orig.ident".
#'
#' @return A ggplot object that visualizes the proportion of each cell type or specified metadata category across the samples or conditions. The plot is a bar plot where the x-axis represents the sample identifiers, the y-axis represents the proportion of each category, and the bars are colored by cell type or metadata category.
#' @export
#'
#' @examples
#' # Assuming 'seurat_obj' is a Seurat object with default metadata columns
#' r2prop_plot(seurat_object = seurat_obj)
#'
#' # Example with specified metadata columns
#' r2prop_plot(seurat_object = seurat_obj, celltype = "cell_type", orig.ident = "sample_name")
#'
#' # Example with different metadata
#' r2prop_plot(seurat_object = seurat_obj, celltype = "cell_type", orig.ident = "condition")
#'
#' @import ggplot2
#' @import dplyr
r2prop_plot <- function(seurat_object, celltype = Idents(seurat_object), orig.ident = "orig.ident") {

  # Extract metadata from the Seurat object
  meta <- seurat_object[[]]

  # Select relevant columns based on the arguments
  meta <- meta[, c(celltype, orig.ident)]

  # Group, summarize, and calculate proportions
  prop_data <- meta %>%
    group_by(.data[[orig.ident]], .data[[celltype]]) %>%
    summarise(count = n(), .groups = 'drop') %>%
    mutate(proportion = count / sum(count) * 100)

  # Plot the data using ggplot2
  ggplot(prop_data, aes(x = .data[[orig.ident]], y = proportion, fill = .data[[celltype]])) +
    geom_bar(stat = "identity") +
    labs(y = "Proportion (%)", x = "", fill = "Celltype")
}
