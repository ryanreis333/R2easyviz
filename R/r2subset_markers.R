#' Subset Marker Genes
#'
#' The `r2subset_markers` function subsets marker genes from a data frame,
#' selecting the top `nfeatures` based on a specified ranking metric.
#'
#' @param markers.result A data frame containing marker gene results (e.g., from `FindAllMarkers`).
#' @param arrange_by A string for the column to order genes by (e.g., `"p_val"`, `"avg_logFC"`).
#'   If `"dif"`, it's calculated as `pct.1 - pct.2`. Default is `"dif"`.
#' @param nfeatures An integer for the number of top genes to select per cluster. Default is 10.
#' @param clusters A string for the column with cluster identifiers. Default is `"cluster"`.
#' @param order A string for sorting order: `"increasing"` or `"decreasing"`. Default is `"decreasing"`.
#'
#' @return A data frame with the top marker genes for each cluster.
#' @export
#'
#' @examples
#' \dontrun{
#' # Create a dummy marker results data frame
#' markers_df <- data.frame(
#'   cluster = rep(c("Cluster1", "Cluster2"), each = 20),
#'   gene = paste0("Gene", 1:40),
#'   avg_logFC = rnorm(40, 0.8, 0.3),
#'   pct.1 = runif(40, 0.4, 0.8),
#'   pct.2 = runif(40, 0.1, 0.5),
#'   p_val_adj = runif(40, 0, 0.05)
#' )
#'
#' # Subset top 5 markers by avg_logFC
#' top_markers <- r2subset_markers(markers.result = markers_df, arrange_by = "avg_logFC", nfeatures = 5)
#'
#' # Subset top 10 markers by 'dif'
#' top_dif_markers <- r2subset_markers(markers.result = markers_df, nfeatures = 10)
#' }
#'
#' @importFrom dplyr group_by arrange slice_head mutate desc
#' @importFrom rlang sym .data
#' @importFrom utils globalVariables
#'
r2subset_markers <- function(markers.result, arrange_by = "dif", nfeatures = 10, clusters = "cluster", order = "decreasing") {
  # Validate order
  if (!order %in% c("increasing", "decreasing")) {
    stop("Parameter 'order' must be either 'increasing' or 'decreasing'.")
  }

  # Check and compute 'dif' if needed
  if (!(arrange_by %in% colnames(markers.result))) {
    if (arrange_by == "dif" && all(c("pct.1", "pct.2") %in% colnames(markers.result))) {
      markers.result <- markers.result %>%
        mutate(dif = .data$pct.1 - .data$pct.2)
    } else {
      stop(paste("Column", arrange_by, "not found in 'markers.result' and cannot be calculated."))
    }
  }

  # Perform grouping and sorting
  markers.subset <- markers.result %>%
    group_by(!!sym(clusters)) %>%
    {
      if (order == "decreasing") {
        arrange(., desc(.data[[arrange_by]]))
      } else {
        arrange(., .data[[arrange_by]])
      }
    } %>%
    slice_head(n = nfeatures)

  return(markers.subset)
}
