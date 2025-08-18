#' Proportional Plot of Cell Types Across Samples
#'
#' This function generates a bar plot showing the proportion of each cell type
#' (or other metadata category) across different samples or conditions from a Seurat object.
#'
#' @param seurat_obj A Seurat object.
#' @param celltype A string specifying the metadata column for cell types.
#'   If `NULL` (default), the active identities (`Idents()`) will be used.
#' @param group.by A string specifying the metadata column for grouping. Default is `"orig.ident"`.
#' @param reorder A logical value. If `TRUE`, samples are reordered based on the
#'   proportion of the most dominant cell type. Default is `FALSE`.
#'
#' @return A ggplot object visualizing the proportions.
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
#' seurat_obj$cell_type <- sample(c("T-cell", "B-cell", "Macrophage"), 20, replace = TRUE)
#' Idents(seurat_obj) <- seurat_obj$cell_type
#'
#' # Plot proportions using active identities
#' r2prop_plot(seurat_obj, group.by = "orig.ident")
#'
#' # Plot proportions using a metadata column
#' r2prop_plot(seurat_obj, celltype = "cell_type", group.by = "orig.ident", reorder = TRUE)
#' }
#'
#' @import ggplot2
#' @importFrom Seurat FetchData Idents
#' @importFrom dplyr %>% group_by summarise mutate pull
#' @importFrom rlang .data
#' @importFrom utils globalVariables
#'
r2prop_plot <- function(seurat_obj, celltype = NULL, group.by = "orig.ident", reorder = FALSE) {

  # If celltype is NULL, use Idents
  if (is.null(celltype)) {
    meta <- data.frame(
      "celltype" = Idents(seurat_obj),
      "group" = FetchData(seurat_obj, vars = group.by)[[1]]
    )
    celltype_col <- "celltype"
    group_col <- "group"
  } else {
    # Check if the specified columns exist in the metadata
    if (!all(c(celltype, group.by) %in% colnames(seurat_obj@meta.data))) {
      stop("One or both specified columns not found in Seurat object metadata.")
    }
    meta <- FetchData(seurat_obj, vars = c(celltype, group.by))
    celltype_col <- celltype
    group_col <- group.by
  }

  # Calculate proportions
  prop_data <- meta %>%
    group_by(.data[[group_col]], .data[[celltype_col]]) %>%
    summarise(count = n(), .groups = 'drop') %>%
    group_by(.data[[group_col]]) %>%
    mutate(proportion = count / sum(count) * 100)

  # Reorder if requested
  if (reorder) {
    dominant_order <- prop_data %>%
      group_by(.data[[group_col]]) %>%
      summarise(max_proportion = max(.data$proportion)) %>%
      arrange(desc(.data$max_proportion)) %>%
      pull(.data[[group_col]])
    prop_data[[group_col]] <- factor(prop_data[[group_col]], levels = dominant_order)
  }

  # Create the plot
  p <- ggplot(prop_data, aes(x = .data[[group_col]], y = .data$proportion, fill = .data[[celltype_col]])) +
    geom_bar(stat = "identity") +
    labs(y = "Proportion (%)", x = "Sample", fill = "Celltype") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  return(p)
}
