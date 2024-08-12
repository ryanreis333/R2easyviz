#' Title
#'
#' @param seurat_object A seurat object you want to create the plot from
#' @param celltype The celltype or other metadata you want to calculate proportions from
#' @param orig.ident The orig.ident or sample name you want to split the data by
#'
#' @return A ggplot that shows the proportion of each celltype across samples
#' @export
#'
#' @examples r2propplot()
#'
r2propplot <- function(seurat_object, celltype, orig.ident, disease) {

  # Extract metadata from the Seurat object
  meta <- seurat_object[[]]

  # Select relevant columns based on the arguments
  meta <- meta[, c(celltype, orig.ident, disease)]

  # Group, summarize, and calculate proportions
  prop_data <- meta %>%
    group_by(.data[[orig.ident]], .data[[celltype]]) %>%
    summarise(count = n()) %>%
    mutate(proportion = count / sum(count) * 100)

  # Plot the data using ggplot2
  ggplot(prop_data, aes(x = .data[[orig.ident]], y = proportion, fill = .data[[celltype]])) +
    geom_bar(stat = "identity") +
    labs(y = "Proportion (%)", x = "", fill = "Cell Type")
}
