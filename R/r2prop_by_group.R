#' Grouped Bar Plot with Individual Sample Points for Seurat Object
#'
#' @param seurat_object A Seurat object containing the data to be plotted.
#' @param celltype A string specifying the column name for cell types.
#' @param group A string specifying the column name for grouping.
#' @param orig.ident A string specifying the column name for sample identifiers.
#'
#' @return A ggplot object that visualizes the grouped bar plot with individual sample points.
#' @export
#'
#' @import ggplot2
#' @import dplyr
r2prop_by_group <- function(seurat_object, celltype = Idents(seurat_object), group = "group", orig.ident = "orig.ident") {
  # Extract metadata from the Seurat object
  meta <- seurat_object[[]]

  # Calculate proportions for each sample
  sample_prop_data <- meta %>%
    group_by(.data[[orig.ident]], .data[[celltype]], .data[[group]]) %>%
    summarise(count = n(), .groups = 'drop') %>%
    group_by(.data[[orig.ident]], .data[[group]]) %>%
    mutate(proportion = count / sum(count) * 100) %>%
    ungroup()

  # Calculate mean proportions for each group and cell type
  summary_data <- sample_prop_data %>%
    group_by(.data[[celltype]], .data[[group]]) %>%
    summarise(mean_proportion = mean(proportion), .groups = 'drop')

  # Create the plot
  ggplot(sample_prop_data, aes(x = .data[[celltype]], y = proportion, fill = .data[[group]], color = .data[[group]])) +
    geom_bar(data = summary_data, aes(y = mean_proportion), stat = "identity",
             position = position_dodge(width = 0.9), width = 0.8, alpha = 0.7) +
    geom_point(position = position_dodge(width = 0.9),
               size = 3) +
    scale_fill_manual(values = c("#44BFC7", "#F8C151")) +  # Adjust colors as needed
    scale_color_manual(values = c("#44BFC7", "#F8C151")) +  # Adjust colors as needed
    labs(x = "Cell Type", y = "Proportion (%)", fill = "Group", color = "Group") +
    theme_minimal() +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
  return(summary_data)
  return(sample_prop_data)
}
