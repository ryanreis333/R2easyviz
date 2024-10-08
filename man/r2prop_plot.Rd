% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/r2prop_plot.R
\name{r2prop_plot}
\alias{r2prop_plot}
\title{Proportional Plot of Cell Types Across Samples}
\usage{
r2prop_plot(
  seurat_object,
  celltype = "celltype",
  orig.ident = "orig.ident",
  reorder = FALSE
)
}
\arguments{
\item{seurat_object}{A Seurat object containing the data to be plotted. The object should include metadata columns corresponding to the cell types and sample identifiers.}

\item{celltype}{A string specifying the column name in the Seurat object's metadata that contains the cell type or other categorical metadata to be plotted. By default, it uses "celltype" which refers to the identities of the cells.}

\item{orig.ident}{A string specifying the column name in the Seurat object's metadata that contains the sample identifiers or conditions by which the data should be split. The default is "orig.ident".}

\item{reorder}{A logical value indicating whether to reorder the samples based on the proportion of the most dominant cell type. If `TRUE`, samples will be reordered; if `FALSE` (default), the original order will be used.}
}
\value{
A ggplot object that visualizes the proportion of each cell type or specified metadata category across the samples or conditions. The plot is a bar plot where the x-axis represents the sample identifiers, the y-axis represents the proportion of each category, and the bars are colored by cell type or metadata category.
}
\description{
This function generates a bar plot showing the proportion of each cell type or specified metadata category across different samples or conditions from a Seurat object. By default, it uses the identities of the cells and the "orig.ident" column for grouping and splitting. Optionally, you can reorder the samples based on the proportion of the most dominant cell type.
}
\examples{
# Assuming 'seurat_obj' is a Seurat object with default metadata columns
r2prop_plot(seurat_object = seurat_obj)

# Example with specified metadata columns
r2prop_plot(seurat_object = seurat_obj, celltype = "cell_type", orig.ident = "sample_name")

# Example with different metadata
r2prop_plot(seurat_object = seurat_obj, celltype = "cell_type", orig.ident = "condition")

# Example with reordered samples
r2prop_plot(seurat_object = seurat_obj, celltype = "cell_type", orig.ident = "sample_name", reorder = TRUE)

}
