#' Remove Specific Features from a Seurat Object
#'
#' The `r2remove_features` function is used to remove specific sets of features (genes) from the `VariableFeatures` of a Seurat object. The function allows users to exclude features based on their category (e.g., TCR genes, IG genes, mitochondrial genes) and supports both human and mouse gene patterns.
#'
#' @param seurat_obj A Seurat object containing single-cell RNA-seq data. This object must have a `VariableFeatures` slot.
#' @param features_list A character vector specifying which categories of features to exclude. The available categories are:
#' \itemize{
#'   \item `"TCR_genes"`: T-cell receptor genes
#'   \item `"IG_genes"`: Immunoglobulin genes
#'   \item `"Mito"`: Mitochondrial genes
#'   \item `"Ribo"`: Ribosomal genes
#'   \item `"Cell.Cycle"`: Cell cycle related genes (requires `cc.genes` object)
#'   \item `"Other"`: Other specified genes
#' }
#' Default is all categories: `c("TCR_genes", "IG_genes", "Mito", "Ribo", "Cell.Cycle", "Other")`.
#' @param species A character string specifying the species of the data. Choices are `"human"` or `"mouse"`. Default is `"human"`.
#'
#' @return A Seurat object with the specified features removed from the `VariableFeatures` slot.
#' @export
#'
#' @examples
#' # Remove IG genes and mitochondrial genes from the VariableFeatures of a Seurat object for human data
#' seurat_filtered <- r2remove_features(seurat_obj = pbmc, features_list = c("IG_genes", "Mito"))
#'
#' # Remove TCR genes and ribosomal genes from the VariableFeatures of a Seurat object for mouse data
#' seurat_filtered <- r2remove_features(seurat_obj = mouse_seurat, features_list = c("TCR_genes", "Ribo"), species = "mouse")
#'
#' @seealso `VariableFeatures` to view or set the variable features in a Seurat object.
#' @importFrom Seurat VariableFeatures
#'
r2remove_features <- function(seurat_obj,
                              features_list = c("TCR_genes", "IG_genes", "Mito", "Ribo", "Cell.Cycle", "Other"),
                              species = "human") {

  # Extract all gene names from the Seurat object
  all.genes <- rownames(seurat_obj)

  # Initialize empty vectors for each category
  IG_genes <- TCR_genes <- Mito <- Ribo <- Cell.Cycle <- Other <- character(0)

  # Human gene patterns
  if (species == "human") {
    if ("IG_genes" %in% features_list) {
      IG_genes <- c(grep("^IGJ", all.genes, value=TRUE),
                    grep("^IGH", all.genes, value=TRUE),
                    grep("^IGK", all.genes, value=TRUE),
                    grep("^IGL", all.genes, value=TRUE))
    }
    if ("TCR_genes" %in% features_list) {
      TCR_genes <- c(grep("^TRA", all.genes, value=TRUE),
                     grep("^TRB", all.genes, value=TRUE),
                     grep("^TRD", all.genes, value=TRUE),
                     grep("^TRG", all.genes, value=TRUE))
    }
    if ("Mito" %in% features_list) {
      Mito <- grep("^MT-", all.genes, value=TRUE)
    }
    if ("Ribo" %in% features_list) {
      Ribo <- c(grep("^RPL", all.genes, value=TRUE),
                grep("^RPS", all.genes, value=TRUE))
    }
    if ("Cell.Cycle" %in% features_list) {
      Cell.Cycle <- c(cc.genes$s.genes, cc.genes$g2m.genes)
    }
    if ("Other" %in% features_list) {
      Other <- c(grep("^MTMR", all.genes, value=TRUE),
                 grep("^MTND", all.genes, value=TRUE),
                 "NEAT1", "TMSB4X", "TMSB10")
    }
  }

  # Mouse gene patterns
  if (species == "mouse") {
    if ("IG_genes" %in% features_list) {
      IG_genes <- c(grep("^Igj", all.genes, value=TRUE),
                    grep("^Igh", all.genes, value=TRUE),
                    grep("^Igk", all.genes, value=TRUE),
                    grep("^Igl", all.genes, value=TRUE))
    }
    if ("TCR_genes" %in% features_list) {
      TCR_genes <- c(grep("^Tra", all.genes, value=TRUE),
                     grep("^Trb", all.genes, value=TRUE),
                     grep("^Trd", all.genes, value=TRUE),
                     grep("^Trg", all.genes, value=TRUE))
    }
    if ("Mito" %in% features_list) {
      Mito <- grep("^mt-", all.genes, value=TRUE)
    }
    if ("Ribo" %in% features_list) {
      Ribo <- c(grep("^Rpl", all.genes, value=TRUE),
                grep("^Rps", all.genes, value=TRUE))
    }
    if ("Cell.Cycle" %in% features_list) {
      Cell.Cycle <- c(cc.genes$s.genes, cc.genes$g2m.genes)
    }
    if ("Other" %in% features_list) {
      Other <- c(grep("^Mtmr", all.genes, value=TRUE),
                 grep("^Mtnd", all.genes, value=TRUE),
                 "Neat1", "Tmsb4x", "Tmsb10")
    }
  }

  # Combine all genes to exclude based on features_list
  exclude_genes <- c(IG_genes, TCR_genes, Mito, Ribo, Cell.Cycle, Other)

  # Remove the specified genes from the VariableFeatures
  VariableFeatures(seurat_obj) <- setdiff(VariableFeatures(seurat_obj), exclude_genes)

  return(seurat_obj)
}
