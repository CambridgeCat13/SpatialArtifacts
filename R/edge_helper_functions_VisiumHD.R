#' Fill isolated spots in `3x3` focal window
#'
#' Fills a center pixel if it is completely surrounded by outlier spots.
#'
#' @param x Numeric vector of length 9 representing a `3x3` window where \code{x[5]} is center
#' @return Integer: 1 if center should be filled, 0 otherwise, NA if center is NA
#' @keywords internal
#' @noRd
my_fill <- function(x) {
  center <- x[5]
  if(is.na(center)) return(NA_integer_)
  if(center == 1) return(1)
  
  # Treat NA as outlier (conservative approach)
  x[is.na(x)] <- 1
  
  # If all 8 neighbors are outliers, fill center
  if(sum(x[-5]) == 8) return(1)
  else return(0)
}


#' Fill spots outlined by outliers in `5x5` window
#'
#' Fills a center pixel if completely outlined by outliers in `5x5` perimeter.
#'
#' @param x Numeric vector of length 25 representing a `5x5` window where \code{x[13]} is center
#' @return Integer: 1 if center should be filled, 0 otherwise, NA if center is NA
#' @keywords internal
#' @noRd
my_outline <- function(x) {
  center <- x[13]
  if(is.na(center)) return(NA_integer_)
  if(center == 1) return(1)
  x[is.na(x)] <- 1
  # Check 16 border pixels: top row (1-5), bottom row (21-25),
  # left column (6,11,16), right column (10,15,20)
  indices <- c(1:5, 21:25, 6, 11, 16, 10, 15, 20)
  sum_edges <- sum(x[indices])
  if(sum_edges == 16) return(1)
  else return(0)
}


#' Fill spots with outliers in cardinal directions (star pattern)
#'
#' Fills center if all four cardinal neighbors (N, S, E, W) are outliers.
#'
#' @param x Numeric vector of length 9 in star pattern arrangement
#' @return Integer: 1 if center should be filled, 0 otherwise, NA if center is NA
#' @keywords internal
#' @noRd
my_fill_star <- function(x) {
  center <- x[5]
  if(is.na(center)) return(NA_integer_)
  if(center == 1) return(1)
  
  # Check cardinal directions: positions 2 (N), 4 (W), 6 (E), 8 (S)
  neighbors <- x[c(2, 4, 6, 8)]
  if(sum(neighbors, na.rm = TRUE) == 4) return(1)
  else return(0)
}


#' Morphological transformations using terra (optimized for VisiumHD)
#'
#' Terra-only implementation of focal operations for connecting outlier regions.
#'
#' @param r Terra SpatRaster with binary values (1=outlier, 0=normal)
#' @param min_cluster_size Minimum size for small hole removal
#'
#' @return Processed SpatRaster
#'
#' @keywords internal
#' @noRd
focal_transformations_terra <- function(r, min_cluster_size = 5) {
  
  if (!inherits(r, "SpatRaster")) {
    stop("This function requires a terra::SpatRaster object.\n",
         "For standard Visium, use the raster-based version.")
  }
  
  # Step 1: 3*3 fill (surrounded by outliers)
  r2 <- terra::focal(r, w = 3, fun = my_fill, na.policy = "all")
  
  # Step 2: 5*5 outline (outlined by outliers)
  r3 <- terra::focal(r2, w = 5, fun = my_outline, na.policy = "all")
  
  # Step 3: Star pattern (cardinal directions)
  star_m <- matrix(c(0, 1, 0, 1, 1, 1, 0, 1, 0), 3, 3)
  r3_s <- terra::focal(r3, w = star_m, fun = my_fill_star, na.policy = "all")
  
  # Step 4: Remove small holes (isolated normal regions)
  # Invert: 1->0, 0->1
  rev_r3 <- terra::ifel(r3_s == 1, 0, 1)
  rev_c3 <- terra::patches(rev_r3, directions = 8, zeroAsNA = TRUE)
  rev_c3_values <- terra::values(rev_c3)
  
  if (!is.null(rev_c3_values) && any(!is.na(rev_c3_values))) {
    # Count cluster sizes
    tbl <- table(rev_c3_values[!is.na(rev_c3_values)])
    
    if (length(tbl) > 0) {
      # Find small clusters (holes) to fill
      flip_clump <- as.numeric(names(tbl)[tbl < min_cluster_size])
      
      if (length(flip_clump) > 0) {
        # Fill small holes by setting them to outlier
        r4 <- r3_s
        r4_values <- terra::values(r4)
        r4_values[rev_c3_values %in% flip_clump] <- 1
        terra::values(r4) <- r4_values
        return(r4)
      }
    }
  }
  
  return(r3_s)
}


#' Efficient coordinate mapping without intermediate raster
#'
#' Maps cluster coordinates back to original spot identifiers using
#' direct matching instead of creating intermediate raster objects.
#'
#' @param xyz_orig Original data frame with x, y coordinates and rownames
#' @param cluster_df Cluster data frame with x, y, cluster_id, size
#'
#' @return Data frame with spotcode, clumpID, clumpSize
#'
#' @keywords internal
#' @noRd
map_coordinates_to_spots_efficient <- function(xyz_orig, cluster_df) {
  # Create coordinate lookup table from original data
  coord_lookup <- data.frame(
    x = xyz_orig[, 1],
    y = xyz_orig[, 2],
    spotcode = rownames(xyz_orig),
    stringsAsFactors = FALSE
  )
  # Round coordinates to avoid floating point issues
  coord_lookup$x <- round(coord_lookup$x, 6)
  coord_lookup$y <- round(coord_lookup$y, 6)
  cluster_df$x <- round(cluster_df$x, 6)
  cluster_df$y <- round(cluster_df$y, 6)
  
  # Merge cluster assignments with spot identifiers
  result <- merge(
    cluster_df,
    coord_lookup,
    by = c("x", "y"),
    all.x = FALSE, 
    all.y = FALSE
  )
  
  if (nrow(result) == 0) {
    return(data.frame(
      spotcode = character(0),
      clumpID = character(0),
      clumpSize = numeric(0)
    ))
  }
  
  output <- data.frame(
    spotcode = result$spotcode,
    clumpID = result$cluster_id,
    clumpSize = result$size,
    stringsAsFactors = FALSE
  )
  
  return(output)
}

#' Detect problem areas in VisiumHD data using morphological processing
#'
#' Terra-optimized implementation for VisiumHD. Uses morphological operations 
#' to connect outlier regions and identify connected components.
#'
#' @param .xyz Data frame with columns: x, y, outlier (binary 0/1)
#' @param uniqueIdentifier Character string for cluster ID prefix
#' @param min_cluster_size Minimum cluster size in bins (default: 5)
#' @param resolution VisiumHD resolution: "8um" or "16um" (default: "8um")
#'
#' @return Data frame with columns: spotcode, clumpID, clumpSize
#' 
#' @examples
#' library(terra)
#' 
#' # 1. Create a mock VisiumHD-like coordinate dataframe
#' # 10x10 grid
#' coords <- expand.grid(x = 1:10, y = 1:10)
#' .xyz <- data.frame(
#'   x = coords$x,
#'   y = coords$y,
#'   outlier = 0
#' )
#' rownames(.xyz) <- paste0("bin_", 1:100)
#' 
#' # 2. Create a "problem area" (cluster of outliers)
#' # Make a 3x3 block of outliers in the center
#' center_mask <- .xyz$x >= 4 & .xyz$x <= 6 & .xyz$y >= 4 & .xyz$y <= 6
#' .xyz$outlier[center_mask] <- 1
#' 
#' # 3. Run detection
#' clusters <- problemAreas_WithMorphology_terra(
#'   .xyz, 
#'   uniqueIdentifier = "TEST", 
#'   min_cluster_size = 2,
#'   resolution = "16um"
#' )
#' 
#' # Check results
#' print(clusters)
#' 
#' @export
problemAreas_WithMorphology_terra <- function(.xyz,
                                              uniqueIdentifier = NA,
                                              min_cluster_size = 5,
                                              resolution = "8um") {
  
  # Validate inputs
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required for VisiumHD processing")
  }
  
  if (ncol(.xyz) < 3) {
    stop("'.xyz' must have at least 3 columns: x, y, outlier")
  }
  
  # Count outliers
  n_outliers <- sum(.xyz[, 3] == 1, na.rm = TRUE)
  
  if (n_outliers == 0) {
    message("    No outliers to cluster")
    return(data.frame(
      spotcode = character(0),
      clumpID = character(0),
      clumpSize = numeric(0)
    ))
  }
  
  message(sprintf("    Clustering %d outlier bins...", n_outliers))
  # STEP 1: CREATE RASTER
  # Get coordinate ranges
  min_x <- min(.xyz[, 1])
  max_x <- max(.xyz[, 1])
  min_y <- min(.xyz[, 2])
  max_y <- max(.xyz[, 2])
  
  # Create raster with exact coordinate system
  nrows <- max_y - min_y + 1
  ncols <- max_x - min_x + 1
  
  message(sprintf("    Creating raster: %d rows * %d cols", nrows, ncols))
  r <- terra::rast(
    nrows = nrows,
    ncols = ncols,
    xmin = min_x - 0.5,  # Cell centers at integer coordinates
    xmax = max_x + 0.5,
    ymin = min_y - 0.5,
    ymax = max_y + 0.5,
    vals = 0
  )
  
  # Set outlier cells to 1
  for (i in 1:nrow(.xyz)) {
    if (.xyz[i, 3] == 1) {
      # rowFromY and colFromX handle coordinate-to-cell mapping
      row_idx <- terra::rowFromY(r, .xyz[i, 2])
      col_idx <- terra::colFromX(r, .xyz[i, 1])
      cell_idx <- terra::cellFromRowCol(r, row_idx, col_idx)
      r[cell_idx] <- 1
    }
  }
  # STEP 2: MORPHOLOGICAL TRANSFORMATIONS
  
  message("    Applying morphological transformations...")
  
  # Adjust min_cluster_size for resolution
  if (resolution == "16um") {
    # 16um has 4* fewer bins than 8um (2* in each dimension)
    # So reduce min_cluster_size accordingly
    adj_min_size <- max(2, min_cluster_size %/% 4)
    message(sprintf("    [16um Mode] Adjusted min_cluster_size: %d -> %d",
                    min_cluster_size, adj_min_size))
    min_cluster_size <- adj_min_size
  }
  
  # Apply sequential morphological operations
  r_clean <- focal_transformations_terra(
    r,
    min_cluster_size = min_cluster_size
  )
  # STEP 3: CLUSTERING WITH CORRECT PARAMETERS
  message("    Clustering outlier regions...")
  c1 <- terra::patches(
    r_clean,
    directions = 8,
    zeroAsNA = TRUE 
  )
  
  # Extract cluster information
  cluster_cells <- which(!is.na(terra::values(c1)))
  
  if (length(cluster_cells) == 0) {
    message("    No clusters found after morphology")
    return(data.frame(
      spotcode = character(0),
      clumpID = character(0),
      clumpSize = numeric(0)
    ))
  }
  
  cluster_ids <- terra::values(c1)[cluster_cells]
  cluster_coords <- terra::xyFromCell(c1, cluster_cells)
  
  message(sprintf("    Found %d clustered cells", length(cluster_cells)))

  # STEP 4: BUILD CLUSTER DATA FRAME 
  # Get unique clusters and their sizes
  cluster_sizes <- table(cluster_ids)
  tot <- length(cluster_sizes)
  
  message(sprintf("    Identified %d unique clusters", tot))
  
  if (is.na(uniqueIdentifier)) uniqueIdentifier <- "X"
  
  # Create results for each cluster
  res_list <- vector("list", tot)
  cluster_names <- names(cluster_sizes)
  
  for (i in 1:tot) {
    cid <- as.numeric(cluster_names[i])
    cluster_mask <- cluster_ids == cid
    
    res_list[[i]] <- data.frame(
      x = cluster_coords[cluster_mask, "x"],
      y = cluster_coords[cluster_mask, "y"],
      cluster_id = paste(uniqueIdentifier, cid, sep = "_"),
      size = sum(cluster_mask),
      stringsAsFactors = FALSE
    )
  }
  
  res_df <- do.call(rbind, res_list)
  # STEP 5: MAP BACK TO SPOT IDENTIFIERS
  pAreas <- map_coordinates_to_spots_efficient(.xyz, res_df)
  
  if (nrow(pAreas) == 0) {
    message("    No spots matched to clusters")
    return(data.frame(
      spotcode = character(0),
      clumpID = character(0),
      clumpSize = numeric(0)
    ))
  }

  # STEP 6: FILTER BY MINIMUM CLUSTER SIZE
  # Count actual spots per cluster
  cluster_sizes_in_pAreas <- table(pAreas$clumpID)
  large_clusters <- names(cluster_sizes_in_pAreas)[
    cluster_sizes_in_pAreas >= min_cluster_size
  ]
  pAreas_filtered <- pAreas[pAreas$clumpID %in% large_clusters, ]
  n_clusters <- length(large_clusters)
  n_bins <- nrow(pAreas_filtered)
  n_filtered <- nrow(pAreas) - n_bins
  message(sprintf("    Kept %d clusters (>=%d bins): %d bins total",
                  n_clusters, min_cluster_size, n_bins))
  
  if (n_filtered > 0) {
    n_small_clusters <- tot - n_clusters
    message(sprintf("    Filtered %d bins in %d small clusters",
                    n_filtered, n_small_clusters))
  }
  
  return(pAreas_filtered)
}