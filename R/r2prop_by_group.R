#' Grouped Bar Plot with Individual Sample Points for Seurat Object
#'
#' This function generates a grouped bar plot showing the proportion of each cell type across different groups, with individual sample points and significance markers. It calculates proportions, performs statistical tests (one-way ANOVA or t-test), and adds significance markers to the plot.
#'
#' @param seurat_object A Seurat object containing the data to be plotted. The object should include metadata columns corresponding to the cell types, grouping, and sample identifiers.
#' @param celltype A string specifying the column name in the Seurat object's metadata that contains the cell type or other categorical metadata to be plotted. By default, it uses `Idents(seurat_object)` which refers to the identities of the cells.
#' @param group A string specifying the column name in the Seurat object's metadata that contains the grouping factor for comparing proportions.
#' @param orig.ident A string specifying the column name in the Seurat object's metadata that contains the sample identifiers or conditions by which the data should be split. The default is "orig.ident".
#' @param pvalue_cutoff A numeric value specifying the cutoff for statistical significance. Default is 0.05.
#'
#' @return A list containing:
#' \item{plot}{A ggplot object that visualizes the grouped bar plot with individual sample points. The x-axis represents cell types, the y-axis represents the proportion of each category, and the bars are colored by group. Individual sample points are overlaid on the bars, and asterisks indicate statistical significance.}
#' \item{stats_results}{A data frame containing the p-values and significance markers for each cell type based on the statistical tests performed.}
#' @export
#'
#' @import ggplot2
#' @import dplyr
#' @import RColorBrewer
#'
#' @examples
#' # Assuming 'seurat_obj' is a Seurat object with default metadata columns
#' result <- r2prop_by_group(seurat_object = seurat_obj)
#' print(result$plot)  # Display the plot
#' print(result$stats_results)  # Display the statistical results
#'
#' # Example with specified metadata columns
#' result <- r2prop_by_group(seurat_object = seurat_obj, celltype = "cell_type", group = "group", orig.ident = "sample_name")
#' print(result$plot)  # Display the plot
#' print(result$stats_results)  # Display the statistical results
#'
#' # Example with different p-value cutoff
#' result <- r2prop_by_group(seurat_object = seurat_obj, pvalue_cutoff = 0.01)
#' print(result$plot)  # Display the plot
#' print(result$stats_results)  # Display the statistical results
r2prop_by_group <- function(seurat_object, celltype = "celltype", group = "group", orig.ident = "orig.ident", pvalue_cutoff = 0.05) {

  # Load necessary libraries
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer) # For color palettes

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

  # Perform ANOVA between groups for each cell type, or t-test if 2 groups
  stats_results <- sample_prop_data %>%
    group_by(.data[[celltype]]) %>%
    summarise(
      p_value = {
        tryCatch({
          num_groups <- length(unique(.data[[group]]))
          if (num_groups > 2) {
            anova_result <- aov(proportion ~ .data[[group]])
            summary_result <- summary(anova_result)
            if (length(summary_result[[1]]$`Pr(>F)`) > 0) {
              summary_result[[1]]$`Pr(>F)`[1]
            } else {
              NA
            }
          } else {
            t_test_result <- t.test(proportion ~ .data[[group]])
            t_test_result$p.value
          }
        }, error = function(e) NA)
      },
      .groups = 'drop'
    ) %>%
    mutate(significance = ifelse(!is.na(p_value) & p_value < pvalue_cutoff, "*", ""))

  # Calculate position for asterisks
  asterisk_data <- summary_data %>%
    inner_join(stats_results, by = celltype) %>%  # Ensure significance column is included
    filter(significance != "") %>%
    mutate(y_position = mean_proportion + 3)  # Position slightly above the top of the bars

  # Define color palettes
  n_colors <- length(unique(meta[[group]]))  # Number of unique groups
  color_palette <- switch(as.character(n_colors),
                          "2" = c("#44BFC7", "#F8C151"),  # Custom colors for 2 groups
                          "3" = brewer.pal(3, "Set1"),  # Palette for 3 groups
                          "4" = brewer.pal(4, "Set1"),  # Palette for 4 groups
                          brewer.pal(min(n_colors, 12), "Set1")  # Palette for more than 4 groups
  )

  # Create the bar plot with points for individual samples
  plot <- ggplot(sample_prop_data, aes(x = .data[[celltype]], y = proportion, fill = .data[[group]])) +
    geom_bar(data = summary_data, aes(y = mean_proportion), stat = "identity",
             position = position_dodge(width = 0.9), width = 0.8, alpha = 0.7) +
    geom_point(aes(color = .data[[group]]), position = position_dodge(width = 0.9), size = 3) +
    geom_text(data = asterisk_data, aes(x = .data[[celltype]], y = y_position, label = significance),
              position = position_dodge(width = 0.9), vjust = -0.5, size = 5, color = "black") +  # Adding asterisks above bars
    scale_fill_manual(values = color_palette) +  # Use color palette for fill
    scale_color_manual(values = color_palette, guide = "none") +  # Use color palette for color
    labs(x = "Cell Type", y = "Proportion (%)", fill = "Group") +
    theme_minimal() +
    theme(legend.position = "right",
          legend.key = element_blank(),         # Remove legend key background
          legend.background = element_blank()   # Remove legend background
    )

  return(list(plot = plot, stats_results = stats_results))
}
