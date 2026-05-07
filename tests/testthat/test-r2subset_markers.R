test_that("r2subset_markers errors on invalid order", {
  df <- data.frame(
    cluster = rep(c("A", "B"), each = 3),
    gene = letters[1:6],
    pct.1 = c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4),
    pct.2 = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
  )
  expect_error(
    r2subset_markers(df, order = "sideways"),
    "must be either 'increasing' or 'decreasing'"
  )
})

test_that("r2subset_markers computes 'dif' from pct.1 and pct.2 when missing", {
  df <- data.frame(
    cluster = rep(c("A", "B"), each = 2),
    gene = letters[1:4],
    pct.1 = c(0.9, 0.8, 0.6, 0.5),
    pct.2 = c(0.1, 0.2, 0.3, 0.4)
  )
  out <- r2subset_markers(df, arrange_by = "dif", nfeatures = 1)
  expect_true("dif" %in% colnames(out))
  expect_equal(nrow(out), 2)
})

test_that("r2subset_markers returns top n per cluster", {
  df <- data.frame(
    cluster = rep(c("A", "B"), each = 4),
    gene = letters[1:8],
    avg_logFC = c(3, 2, 1, 0.5, 4, 3, 2, 1)
  )
  out <- r2subset_markers(df, arrange_by = "avg_logFC", nfeatures = 2,
                          clusters = "cluster")
  expect_equal(nrow(out), 4)
})
