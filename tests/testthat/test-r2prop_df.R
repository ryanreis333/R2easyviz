test_that("r2prop_df rejects non-Seurat input", {
  expect_error(
    r2prop_df(seurat_obj = data.frame(), sample_col = "a", celltype_col = "b"),
    "must be a Seurat object"
  )
})

test_that("r2prop_df errors when sample/celltype columns are missing", {
  skip_if_not_installed("SeuratObject")

  counts <- matrix(rpois(20, lambda = 3), nrow = 4,
                   dimnames = list(paste0("g", 1:4), paste0("c", 1:5)))
  seu <- SeuratObject::CreateSeuratObject(counts = counts)

  expect_error(
    r2prop_df(seu, sample_col = "missing_col", celltype_col = "orig.ident"),
    "not found in Seurat metadata"
  )
  expect_error(
    r2prop_df(seu, sample_col = "orig.ident", celltype_col = "missing_col"),
    "not found in Seurat metadata"
  )
})

test_that("r2prop_df returns proportions that sum to 1 within each sample", {
  skip_if_not_installed("SeuratObject")

  counts <- matrix(rpois(40, lambda = 3), nrow = 4,
                   dimnames = list(paste0("g", 1:4), paste0("c", 1:10)))
  seu <- SeuratObject::CreateSeuratObject(counts = counts)
  seu$sample <- rep(c("s1", "s2"), each = 5)
  seu$celltype <- rep(c("A", "B"), times = 5)

  out <- r2prop_df(seu, sample_col = "sample", celltype_col = "celltype")

  by_sample <- tapply(out$proportion, out$sample, sum)
  expect_equal(unname(by_sample), rep(1, length(by_sample)))
})
