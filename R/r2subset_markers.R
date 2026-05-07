#' Subset Marker Genes
#'
#' The `r2subset_markers` function subsets marker genes by selecting the
#' top `nfeatures` based on a specified ranking metric. The function groups
#' marker genes by cluster and orders them within each cluster according to
#' the chosen metric (e.g., differential expression). If the specified
#' metric is `"dif"` and the necessary columns (`pct.1` and `pct.2`) are
#' present, it calculates `dif` as the difference between `pct.1` and
#' `pct.2`. The function returns a subset of the top marker genes for each
#' cluster.
#'
#' @param markers.result A data frame or tibble containing marker gene
#'   results. This should include columns for clusters and the metric to
#'   arrange by.
#' @param arrange_by A character string specifying the column name in
#'   `markers.result` that contains the metric for ordering the marker
#'   genes (e.g., "p_val", "avg_logFC", "dif"). Default is `"dif"`.
#' @param nfeatures An integer specifying the number of top marker genes to
#'   select per cluster. Default is 10.
#' @param clusters A character string specifying the column name in
#'   `markers.result` that contains the cluster identifiers. Default is
#'   `"cluster"`.
#' @param order A character string indicating whether to sort in
#'   `"increasing"` or `"decreasing"` order. Default is `"decreasing"`.
#'
#' @return A data frame or tibble with the top marker genes for each
#'   cluster based on the specified metric.
#' @export
#'
#' @examples
#' \dontrun{
#' # Top 10 marker genes by avg_logFC
#' top_markers <- r2subset_markers(
#'   markers.result = markers_df,
#'   arrange_by = "avg_logFC",
#'   nfeatures = 10
#' )
#'
#' # Top 5 marker genes by p-value (increasing order)
#' top_markers <- r2subset_markers(
#'   markers.result = markers_df,
#'   arrange_by = "p_val",
#'   nfeatures = 5,
#'   order = "increasing"
#' )
#' }
#'
#' @importFrom dplyr group_by arrange slice_head mutate desc
#' @importFrom rlang sym expr
r2subset_markers <- function(markers.result, arrange_by = "dif", nfeatures = 10, clusters = "cluster", order = "decreasing") {
  # Validate order
  if (!order %in% c("increasing", "decreasing")) {
    stop("Parameter 'order' must be either 'increasing' or 'decreasing'.")
  }

  # Check and compute 'dif' if needed
  if (!(arrange_by %in% colnames(markers.result))) {
    if (arrange_by == "dif" && all(c("pct.1", "pct.2") %in% colnames(markers.result))) {
      markers.result <- markers.result %>%
        dplyr::mutate(dif = .data$pct.1 - .data$pct.2)
    } else {
      stop(paste("Column", arrange_by, "not found in 'markers.result'."))
    }
  }

  # Build the arrange expression (descending or ascending) without bare `.`
  arrange_expr <- if (order == "decreasing") {
    rlang::expr(dplyr::desc(!!rlang::sym(arrange_by)))
  } else {
    rlang::sym(arrange_by)
  }

  # Perform grouping and sorting
  markers.subset <- markers.result %>%
    dplyr::group_by(!!rlang::sym(clusters)) %>%
    dplyr::arrange(!!arrange_expr) %>%
    dplyr::slice_head(n = nfeatures)

  return(markers.subset)
}
