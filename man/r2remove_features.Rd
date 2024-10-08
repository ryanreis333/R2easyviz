% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/r2remove_features.R
\name{r2remove_features}
\alias{r2remove_features}
\title{Remove Specific Features from a Seurat Object}
\usage{
r2remove_features(
  seurat_obj,
  features_list = c("TCR_genes", "IG_genes", "Mito", "Ribo", "Cell.Cycle", "Other"),
  species = "human"
)
}
\arguments{
\item{seurat_obj}{A Seurat object containing single-cell RNA-seq data. This object must have a `VariableFeatures` slot.}

\item{features_list}{A character vector specifying which categories of features to exclude. The available categories are:
\itemize{
  \item `"TCR_genes"`: T-cell receptor genes
  \item `"IG_genes"`: Immunoglobulin genes
  \item `"Mito"`: Mitochondrial genes
  \item `"Ribo"`: Ribosomal genes
  \item `"Cell.Cycle"`: Cell cycle related genes (requires `cc.genes` object)
  \item `"Other"`: Other specified genes
}
Default is all categories: `c("TCR_genes", "IG_genes", "Mito", "Ribo", "Cell.Cycle", "Other")`.}

\item{species}{A character string specifying the species of the data. Choices are `"human"` or `"mouse"`. Default is `"human"`.}
}
\value{
A Seurat object with the specified features removed from the `VariableFeatures` slot.
}
\description{
The `r2remove_features` function is used to remove specific sets of features (genes) from the `VariableFeatures` of a Seurat object. The function allows users to exclude features based on their category (e.g., TCR genes, IG genes, mitochondrial genes) and supports both human and mouse gene patterns.
}
\examples{
# Remove IG genes and mitochondrial genes from the VariableFeatures of a Seurat object for human data
seurat_filtered <- r2remove_features(seurat_obj = pbmc, features_list = c("IG_genes", "Mito"))

# Remove TCR genes and ribosomal genes from the VariableFeatures of a Seurat object for mouse data
seurat_filtered <- r2remove_features(seurat_obj = mouse_seurat, features_list = c("TCR_genes", "Ribo"), species = "mouse")

}
\seealso{
`VariableFeatures` to view or set the variable features in a Seurat object.
}
