#' Subset Marker Genes
#'
#' The `r2subset_markers` function subsets marker genes by selecting the top `nfeatures` based on a specified ranking metric. The function groups marker genes by cluster and orders them within each cluster according to the chosen metric (e.g., differential expression). If the specified metric is `"dif"` and the necessary columns (`pct.1` and `pct.2`) are present, it calculates `dif` as the difference between `pct.1` and `pct.2`. The function returns a subset of the top marker genes for each cluster.
#'
#' @param markers.result A data frame or tibble containing marker gene results. This should include columns for clusters and the metric to arrange by.
#' @param arrange_by A character string specifying the column name in `markers.result` that contains the metric for ordering the marker genes (e.g., "p_val", "avg_logFC", "dif"). Default is `"dif"`.
#' @param nfeatures An integer specifying the number of top marker genes to select per cluster. Default is 10.
#' @param clusters A character string specifying the column name in `markers.result` that contains the cluster identifiers. Default is `"cluster"`.
#'
#' @return A data frame or tibble with the top marker genes for each cluster based on the specified metric.
#' @export
#'
#' @examples
#' # Subset the top 10 marker genes by differential expression from a result data frame
#' top_markers <- r2subset_markers(markers.result = markers_df, arrange_by = "avg_logFC", nfeatures = 10)
#'
#' # Subset the top 5 marker genes based on p-value
#' top_markers <- r2subset_markers(markers.result = markers_df, arrange_by = "p_val", nfeatures = 5)
#'
#' # Subset the top 10 marker genes based on a custom 'dif' metric, calculated as pct.1 - pct.2
#' top_markers <- r2subset_markers(markers.result = markers_df, arrange_by = "dif", nfeatures = 10)
#'
#' @importFrom dplyr group_by arrange slice_head mutate
#' @importFrom rlang sym
#'
r2subset_markers <- function(markers.result, arrange_by = "dif", nfeatures = 10, clusters = "cluster") {
  # Check if the arrange_by column exists
  if (!(arrange_by %in% colnames(markers.result))) {
    # If arrange_by is "dif" and dif column doesn't exist, calculate it
    if (arrange_by == "dif" && all(c("pct.1", "pct.2") %in% colnames(markers.result))) {
      markers.result <- markers.result %>%
        mutate(dif = pct.1 - pct.2)
      arrange_by <- "dif"  # Ensure we're using the newly created 'dif' column
    } else {
      stop(paste("Column", arrange_by, "not found in 'markers.result'."))
    }
  }

  # Subset marker genes by grouping and arranging
  markers.subset <- markers.result %>%
    group_by(!!sym(clusters)) %>%
    arrange(desc(!!sym(arrange_by))) %>%
    slice_head(n = nfeatures)

  return(markers.subset)
}

