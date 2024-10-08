% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/r2heatmap.R
\name{r2heatmap}
\alias{r2heatmap}
\title{Generate a Heatmap of Marker Genes}
\usage{
r2heatmap(
  seurat_obj,
  FindAllMarkersObj,
  group_by = "celltype",
  ncells = 500,
  viridis_color = TRUE,
  nfeatures = 5,
  arrange_by = "dif",
  barcode_column = "barcodes"
)
}
\arguments{
\item{seurat_obj}{A `Seurat` object containing the single-cell RNA-seq data. This object should include metadata with cell grouping information.}

\item{FindAllMarkersObj}{A data frame or tibble containing the results from `FindAllMarkers`, which includes information about marker genes. This should have columns for gene names, clusters, and the metric used for ranking.}

\item{group_by}{A character string specifying the metadata column used for grouping cells in the heatmap. Default is `"celltype"`.}

\item{ncells}{An integer specifying the number of cells to sample per group for the heatmap. If `ncells` is greater than the minimum number of cells in any group, it will be adjusted to the minimum. Default is `500`.}

\item{viridis_color}{A logical value indicating whether to use the Viridis color scale for the heatmap. If `TRUE`, the Viridis color scale will be applied; otherwise, a default color scale will be used. Default is `TRUE`.}

\item{nfeatures}{An integer specifying the number of top marker genes to include in the heatmap per cluster. Default is `5`.}

\item{arrange_by}{A character string specifying the column name in `FindAllMarkersObj` used for ordering marker genes. If `"dif"`, and the `dif` column is not present, it will be calculated as the difference between `pct.1` and `pct.2`. Default is `"dif"`.}

\item{barcode_column}{A character string specifying the column name in the metadata that contains barcode information. Default is `"barcodes"`.}
}
\value{
A `ggplot` object representing the heatmap of the selected marker genes.
}
\description{
The `r2heatmap` function creates a heatmap of the top marker genes for a given `Seurat` object. It allows you to visualize the expression of selected marker genes across cells, with options for color scaling and grouping. The function includes checks for column existence and dynamically adjusts the number of cells sampled based on the available data.
}
\examples{
# Generate a heatmap for the top 5 marker genes, using the Viridis color scale
heatmap_plot <- r2heatmap(seurat_obj = pbmc, FindAllMarkersObj = markers_df, group_by = "celltype", nfeatures = 5)

# Generate a heatmap for the top 10 marker genes, without using the Viridis color scale
heatmap_plot <- r2heatmap(seurat_obj = pbmc, FindAllMarkersObj = markers_df, group_by = "celltype", nfeatures = 10, viridis_color = FALSE)

}
