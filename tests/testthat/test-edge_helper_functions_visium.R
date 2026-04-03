# test-edge_helper_functions_visium.R

#   Tests for my_fill  
test_that("my_fill returns 1 when center is already 1", {
  x <- rep(0, 9)
  x[5] <- 1
  expect_equal(SpatialArtifacts:::my_fill(x), 1)
})

test_that("my_fill returns 1 when all neighbors are 1", {
  x <- rep(1, 9)
  x[5] <- 0
  expect_equal(SpatialArtifacts:::my_fill(x), 1)
})

test_that("my_fill returns 0 when not surrounded", {
  x <- rep(0, 9)
  expect_equal(SpatialArtifacts:::my_fill(x), 0)
})

test_that("my_fill returns NA when center is NA", {
  x <- rep(0, 9)
  x[5] <- NA
  expect_equal(SpatialArtifacts:::my_fill(x), NA_integer_)
})

#   Tests for my_fill_star  
test_that("my_fill_star returns 1 when center is already 1", {
  x <- rep(0, 9)
  x[5] <- 1
  expect_equal(SpatialArtifacts:::my_fill_star(x), 1)
})

test_that("my_fill_star returns 1 when all 4 cardinal neighbors are 1", {
  x <- rep(0, 9)
  x[c(2, 4, 6, 8)] <- 1  # N, W, E, S
  expect_equal(SpatialArtifacts:::my_fill_star(x), 1)
})

test_that("my_fill_star returns NA when center is NA", {
  x <- rep(0, 9)
  x[5] <- NA
  expect_equal(SpatialArtifacts:::my_fill_star(x), NA_integer_)
})

#   Tests for my_outline  
test_that("my_outline returns 1 when center is already 1", {
  x <- rep(0, 25)
  x[13] <- 1
  expect_equal(SpatialArtifacts:::my_outline(x), 1)
})

test_that("my_outline returns 1 when all border pixels are 1", {
  x <- rep(1, 25)
  x[13] <- 0
  expect_equal(SpatialArtifacts:::my_outline(x), 1)
})

test_that("my_outline returns NA when center is NA", {
  x <- rep(1, 25)
  x[13] <- NA
  expect_equal(SpatialArtifacts:::my_outline(x), NA_integer_)
})

#   Tests for focal_transformations  
test_that("focal_transformations returns a SpatRaster", {
  library(terra)
  m <- matrix(0, nrow = 5, ncol = 5)
  m[1, ] <- 1  # top row as outliers
  r <- terra::rast(m)
  result <- focal_transformations(r, min_cluster_size = 3)
  expect_true(inherits(result, "SpatRaster"))
})

test_that("focal_transformations preserves outlier values", {
  library(terra)
  m <- matrix(0, nrow = 5, ncol = 5)
  m[1, ] <- 1
  r <- terra::rast(m)
  result <- focal_transformations(r, min_cluster_size = 3)
  vals <- terra::values(result)
  expect_true(any(vals == 1, na.rm = TRUE))
})

#   Tests for clumpEdges  
test_that("clumpEdges detects edge artifacts", {
  spot_data <- data.frame(
    array_row = rep(1:5, each = 5),
    array_col = rep(1:5, times = 5),
    outlier = FALSE
  )
  rownames(spot_data) <- paste0("spot", seq_len(25))
  spot_data$outlier[spot_data$array_row == 1] <- TRUE

  result <- clumpEdges(
    spot_data,
    offTissue = character(0),
    edge_threshold = 0.5,
    min_cluster_size = 3
  )
  expect_true(length(result) > 0)
  expect_true(all(result %in% rownames(spot_data)))
})

test_that("clumpEdges returns empty when no outliers", {
  spot_data <- data.frame(
    array_row = rep(1:5, each = 5),
    array_col = rep(1:5, times = 5),
    outlier = FALSE
  )
  rownames(spot_data) <- paste0("spot", seq_len(25))

  result <- clumpEdges(
    spot_data,
    offTissue = character(0),
    edge_threshold = 0.5,
    min_cluster_size = 3
  )
  expect_equal(length(result), 0)
})

#   Tests for problemAreas  
test_that("problemAreas detects interior artifacts", {
  spot_data <- data.frame(
    array_row = rep(1:5, each = 5),
    array_col = rep(1:5, times = 5),
    outlier = FALSE
  )
  rownames(spot_data) <- paste0("spot", seq_len(25))
  spot_data$outlier[spot_data$array_row %in% 2:3 &
                      spot_data$array_col %in% 2:3] <- TRUE

  result <- problemAreas(
    spot_data,
    offTissue = character(0),
    uniqueIdentifier = "Sample1",
    min_cluster_size = 3
  )
  expect_true(nrow(result) > 0)
  expect_true(all(c("spotcode", "clumpID", "clumpSize") %in% colnames(result)))
})

test_that("problemAreas returns empty df when no outliers", {
  spot_data <- data.frame(
    array_row = rep(1:5, each = 5),
    array_col = rep(1:5, times = 5),
    outlier = FALSE
  )
  rownames(spot_data) <- paste0("spot", seq_len(25))

  result <- problemAreas(
    spot_data,
    offTissue = character(0),
    uniqueIdentifier = "Sample1",
    min_cluster_size = 3
  )
  expect_equal(nrow(result), 0)
})