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
r2propplot <- function(seurat_object = "so", celltype = "celltype", orig.ident = "orig.ident") {

  meta <- seurat_object[[]]
  meta <- meta[,c(celltype,orig.ident)]

  prop_data <- meta %>%
    group_by(orig.ident, celltype) %>%
    summarise(count = n()) %>%
    mutate(proportion = count / sum(count) * 100)

  ggplot(prop_data, aes(x = orig.ident, y = proportion, fill = celltype)) +
    geom_bar(stat = "identity") +
    labs(y = "Proportion (%)", x = "", fill = "Cell Type") +
    ggtitle("") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 20, face = "bold"),      # Title
      axis.title.x = element_text(size = 16),                  # X-axis label
      axis.title.y = element_text(size = 20),                  # Y-axis label
      axis.text.x = element_text(size = 18),                   # X-axis text
      axis.text.y = element_text(size = 18),                   # Y-axis text
      legend.title = element_text(size = 16),                  # Legend title
      legend.text = element_text(size = 14)                    # Legend text
    ) + NoLegend()
}
