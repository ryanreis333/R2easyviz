#' Proportional Plot of Cell Types Across Samples
#'
#' This function generates a bar plot showing the proportion of each cell
#' type or specified metadata category across different samples or
#' conditions from a Seurat object. By default, it uses the identities of
#' the cells and the "orig.ident" column for grouping and splitting.
#' Optionally, you can reorder the samples based on the proportion of the
#' most dominant cell type.
#'
#' @param seurat_object A Seurat object containing the data to be plotted.
#'   The object should include metadata columns corresponding to the cell
#'   types and sample identifiers.
#' @param celltype A string specifying the column name in the Seurat
#'   object's metadata that contains the cell type or other categorical
#'   metadata to be plotted. By default, it uses "celltype".
#' @param group.by A string specifying the column name in the Seurat
#'   object's metadata that contains the sample identifiers or conditions
#'   by which the data should be split. The default is "orig.ident".
#' @param reorder A logical value indicating whether to reorder the samples
#'   based on the proportion of the most dominant cell type. If `TRUE`,
#'   samples will be reordered; if `FALSE` (default), the original order
#'   will be used.
#'
#' @return A ggplot object that visualizes the proportion of each cell type
#'   or specified metadata category across the samples or conditions.
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming 'seurat_obj' is a Seurat object with default metadata columns
#' r2prop_plot(seurat_object = seurat_obj)
#'
#' # Example with specified metadata columns
#' r2prop_plot(
#'   seurat_object = seurat_obj,
#'   celltype = "cell_type",
#'   group.by = "sample_name"
#' )
#'
#' # Example with reordered samples
#' r2prop_plot(
#'   seurat_object = seurat_obj,
#'   celltype = "cell_type",
#'   group.by = "sample_name",
#'   reorder = TRUE
#' )
#' }
#'
#' @importFrom dplyr group_by summarise mutate arrange desc pull n
#' @importFrom ggplot2 ggplot aes geom_bar labs theme_minimal theme element_text
r2prop_plot <- function(seurat_object, celltype = "celltype", group.by = "orig.ident", reorder = FALSE) {

  # Extract metadata from the Seurat object
  meta <- seurat_object[[]]

  # Select relevant columns based on the arguments
  meta <- meta[, c(celltype, group.by)]

  # Group, summarize, and calculate proportions within each sample (group.by)
  prop_data <- meta %>%
    dplyr::group_by(.data[[group.by]], .data[[celltype]]) %>%
    dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
    dplyr::group_by(.data[[group.by]]) %>%
    dplyr::mutate(proportion = .data$count / sum(.data$count) * 100)

  # Conditionally reorder group.by based on the proportion of the most dominant celltype
  if (reorder) {
    dominant_order <- prop_data %>%
      dplyr::group_by(.data[[group.by]]) %>%
      dplyr::summarise(max_proportion = max(.data$proportion)) %>%
      dplyr::arrange(dplyr::desc(.data$max_proportion)) %>%
      dplyr::pull(.data[[group.by]])

    plot <- ggplot2::ggplot(
      prop_data,
      ggplot2::aes(x = factor(.data[[group.by]], levels = dominant_order),
                   y = .data$proportion,
                   fill = .data[[celltype]])
    ) +
      ggplot2::geom_bar(stat = "identity")
  } else {
    plot <- ggplot2::ggplot(
      prop_data,
      ggplot2::aes(x = .data[[group.by]],
                   y = .data$proportion,
                   fill = .data[[celltype]])
    ) +
      ggplot2::geom_bar(stat = "identity")
  }

  plot +
    ggplot2::labs(y = "Proportion (%)", x = "Sample", fill = "Celltype") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}
