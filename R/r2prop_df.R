library(dplyr)
library(tidyr)

#' Calculate Cell Type Proportions by Sample
#'
#' @param data A data.frame or tibble (e.g., Seurat_obj@meta.data).
#' @param sample_col String. Column name for the sample IDs (e.g., "orig.ident", "Sample").
#' @param celltype_col String. Column name for cell types/clusters (e.g., "seurat_clusters", "celltype").
#' @param group_col String (Optional). Metadata to keep for comparison (e.g., "Condition", "Response").
#'
#' @return A dataframe in long format with counts and proportions.
calculate_cell_proportions <- function(data, sample_col, celltype_col, group_col = NULL) {

  # 1. Calculate basic counts per sample and cell type
  # We use .data[[string]] to handle string variable names in dplyr
  counts_df <- data %>%
    group_by(.data[[sample_col]], .data[[celltype_col]]) %>%
    summarise(n = n(), .groups = 'drop')

  # 2. Fill missing combinations with 0
  # If Sample A has no B-cells, we need a row saying "Sample A, B-cell, n=0"
  # complete() expands the dataframe to include all combinations of sample and cell type
  full_df <- counts_df %>%
    complete(.data[[sample_col]], .data[[celltype_col]], fill = list(n = 0))

  # 3. Calculate Proportions
  prop_df <- full_df %>%
    group_by(.data[[sample_col]]) %>%
    mutate(
      total_cells = sum(n),
      proportion = n / total_cells
    ) %>%
    ungroup()

  # 4. Re-attach Group/Condition Metadata
  if (!is.null(group_col)) {
    # Create a unique mapping of Sample -> Group
    # We assume one sample belongs to only one group
    meta_map <- data %>%
      select(all_of(c(sample_col, group_col))) %>%
      distinct()

    # Join the group info back to the proportion table
    prop_df <- prop_df %>%
      left_join(meta_map, by = sample_col)
  }

  return(prop_df)
}

# ==========================================
# EXAMPLE USAGE
# ==========================================

# 1. Mock Data Generation
# Imagine this is your Seurat@meta.data
mock_meta <- data.frame(
  sample_id = c(rep("S1", 10), rep("S2", 10), rep("S3", 10)),
  cell_type = c(
    rep("T-cell", 8), rep("B-cell", 2),   # S1: Mostly T-cells
    rep("B-cell", 10),                    # S2: All B-cells (T-cell prop should be 0)
    rep("T-cell", 5), rep("B-cell", 5)    # S3: 50/50 split
  ),
  condition = c(rep("Control", 10), rep("Control", 10), rep("Disease", 10))
)

# 2. Run Function
results <- calculate_cell_proportions(
  data = mock_meta,
  sample_col = "sample_id",
  celltype_col = "cell_type",
  group_col = "condition"
)

# 3. View Results
print("--- Calculated Proportions ---")
print(results)

# Example Plotting Code (ggplot2)
# library(ggplot2)
# ggplot(results, aes(x = condition, y = proportion, fill = cell_type)) +
#   geom_boxplot() +
#   facet_wrap(~cell_type, scales = "free") +
#   theme_bw()
