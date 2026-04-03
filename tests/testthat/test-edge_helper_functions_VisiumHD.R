# test-edge_helper_functions_VisiumHD.R

# Tests for focal_transformations_terra 
test_that("focal_transformations_terra returns a SpatRaster", {
  skip_if_not_installed("terra")
  m <- matrix(0, nrow = 5, ncol = 5)
  m[1, ] <- 1
  r <- terra::rast(m)
  result <- SpatialArtifacts:::focal_transformations_terra(r, min_cluster_size = 3)
  expect_true(inherits(result, "SpatRaster"))
})

test_that("focal_transformations_terra preserves outlier values", {
  skip_if_not_installed("terra")
  m <- matrix(0, nrow = 5, ncol = 5)
  m[1, ] <- 1
  r <- terra::rast(m)
  result <- SpatialArtifacts:::focal_transformations_terra(r, min_cluster_size = 3)
  vals <- terra::values(result)
  expect_true(any(vals == 1, na.rm = TRUE))
})

test_that("focal_transformations_terra errors on non-SpatRaster input", {
  m <- matrix(0, nrow = 5, ncol = 5)
  expect_error(
    SpatialArtifacts:::focal_transformations_terra(m, min_cluster_size = 3),
    "SpatRaster"
  )
})

# --- Tests for map_coordinates_to_spots_efficient ---
test_that("map_coordinates_to_spots_efficient returns correct columns", {
  xyz_orig <- data.frame(x = 1:5, y = 1:5)
  rownames(xyz_orig) <- paste0("spot_", 1:5)
  
  cluster_df <- data.frame(
    x = c(1, 2),
    y = c(1, 2),
    cluster_id = c("S1_1", "S1_1"),
    size = c(2, 2),
    stringsAsFactors = FALSE
  )
  
  result <- SpatialArtifacts:::map_coordinates_to_spots_efficient(
    xyz_orig, cluster_df
  )
  
  expect_true(all(c("spotcode", "clumpID", "clumpSize") %in% colnames(result)))
  expect_true(nrow(result) > 0)
})

test_that("map_coordinates_to_spots_efficient returns empty df when no match", {
  xyz_orig <- data.frame(x = 1:5, y = 1:5)
  rownames(xyz_orig) <- paste0("spot_", 1:5)
  
  cluster_df <- data.frame(
    x = c(99, 100),
    y = c(99, 100),
    cluster_id = c("S1_1", "S1_1"),
    size = c(2, 2),
    stringsAsFactors = FALSE
  )
  
  result <- SpatialArtifacts:::map_coordinates_to_spots_efficient(
    xyz_orig, cluster_df
  )
  
  expect_equal(nrow(result), 0)
})

# --- Tests for problemAreas_WithMorphology_terra ---
test_that("problemAreas_WithMorphology_terra returns expected columns", {
  skip_if_not_installed("terra")
  
  coords <- expand.grid(x = 1:10, y = 1:10)
  xyz <- data.frame(x = coords$x, y = coords$y, outlier = 0)
  rownames(xyz) <- paste0("bin_", seq_len(nrow(xyz)))
  
  # Interior artifact
  center_mask <- xyz$x >= 4 & xyz$x <= 6 & xyz$y >= 4 & xyz$y <= 6
  xyz$outlier[center_mask] <- 1
  
  result <- problemAreas_WithMorphology_terra(
    xyz,
    uniqueIdentifier = "TEST",
    min_cluster_size = 2,
    resolution = "16um"
  )
  
  expect_true(all(c("spotcode", "clumpID", "clumpSize") %in% colnames(result)))
})

test_that("problemAreas_WithMorphology_terra returns empty df when no outliers", {
  skip_if_not_installed("terra")
  
  coords <- expand.grid(x = 1:10, y = 1:10)
  xyz <- data.frame(x = coords$x, y = coords$y, outlier = 0)
  rownames(xyz) <- paste0("bin_", seq_len(nrow(xyz)))
  
  result <- problemAreas_WithMorphology_terra(
    xyz,
    uniqueIdentifier = "TEST",
    min_cluster_size = 2,
    resolution = "16um"
  )
  
  expect_equal(nrow(result), 0)
})

test_that("problemAreas_WithMorphology_terra detects interior clusters", {
  skip_if_not_installed("terra")
  
  coords <- expand.grid(x = 1:10, y = 1:10)
  xyz <- data.frame(x = coords$x, y = coords$y, outlier = 0)
  rownames(xyz) <- paste0("bin_", seq_len(nrow(xyz)))
  
  center_mask <- xyz$x >= 4 & xyz$x <= 6 & xyz$y >= 4 & xyz$y <= 6
  xyz$outlier[center_mask] <- 1
  
  result <- problemAreas_WithMorphology_terra(
    xyz,
    uniqueIdentifier = "TEST",
    min_cluster_size = 2,
    resolution = "16um"
  )
  
  expect_true(nrow(result) > 0)
  expect_true(all(result$spotcode %in% rownames(xyz)))
})

test_that("problemAreas_WithMorphology_terra errors with insufficient columns", {
  xyz_bad <- data.frame(x = 1:5, y = 1:5)
  
  expect_error(
    problemAreas_WithMorphology_terra(xyz_bad),
    "3 columns"
  )
})