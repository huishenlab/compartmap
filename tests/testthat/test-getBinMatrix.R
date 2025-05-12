test_that("getBinMatrix", {
  mat <- matrix(5, nrow = 10, ncol = 10)
  mat.na <- mat
  mat.na[1, 1] <- NA

  chrs <- c("chr2", "chr2", "chr1", "chr3")
  gr <- GRanges(Rle(chrs, c(1, 1, 1, 2)), IRanges(1:5, width = 1))
  gr.empty <- GRanges()

  expect_error(
    getBinMatrix(mat.na, gr),
    "Matrix must not contain NAs",
    fixed = TRUE
  )

  expect_error(
    getBinMatrix(mat, gr.empty),
    "Provided GRanges must have length equal to the matrix number of rows",
    fixed = TRUE
  )
})

test_that(".makeBins", {
  chrs <- c("chr2", "chr2", "chr1", "chr3")
  gr <- GRanges(Rle(chrs, c(1, 1, 1, 2)), IRanges(1:5, width = 1))
  compartmap:::.makeBins(1, 100, 10, "chr1")

  bins <- list(
    c(start = 1, end = 100, res = 10, chr = "chr1"),
    c(start = 2, end = 256, res = 6, chr = "chr2"),
    c(start = 5, end = 500, res = 2, chr = "chrX")
  )

  lapply(bins, function(input) {
    start <- as.numeric(input['start'])
    end <- as.numeric(input['end'])
    res <- as.numeric(input['res'])
    chr <- input['chr']
    ex.start <- seq(start, end, by = res)
    ex.end <- c(seq((start + res) - 1, end - 1, by = res), end - 1)
    ex.chr <- rep(chr, length(ex.start)) |> as.factor() |> unname()
    gdf <- as.data.frame(compartmap:::.makeBins(start, end, res, chr))

    expect_equal(gdf$start, ex.start)
    expect_equal(gdf$end, ex.end)
    expect_equal(gdf$seqnames, ex.chr)
  })
})

test_that(".summarizeBins", {
  expected.base <- tapply(iris$Sepal.Length, iris$Species, mean) |> as.numeric()

  get_indices <- function(indices) {
    c(rep(indices[1], 50), rep(indices[2], 50), rep(indices[3], 50))
  }

  tester <- function(test.input) {
    ind <- unlist(test.input$indices)
    ind.full <- get_indices(ind)
    test <- compartmap:::.summarizeBins(
      iris$Sepal.Length,
      test.input$binlength,
      ind.full,
      mean
    )
    expected <- rep(0, test.input$binlength)
    expected[ind] <- expected.base
    expect_equal(test, expected)
  }

  index_set <- list(
    c(1, 2, 3),
    c(1, 5, 10),
    c(5, 50, 500)
  )
  binlength_set <- c(3, 5, 10, 50, 100, 500)
  test_set <- expand.grid(binlength = binlength_set, indices = index_set)
  apply(test_set, 1, tester)
})
