#' Detect Edge Artifacts in VisiumHD Data
#'
#' @param spe SpatialExperiment object
#' @param resolution Resolution: "8um" or "16um" (REQUIRED)
#' @param qc_metric QC metric column (default: "sum_gene")
#' @param samples Sample ID column (default: "sample_id")
#' @param mad_threshold MAD threshold (default: 3)
#' @param buffer_width_um Buffer zone width in micrometers (default: 80)
#'   Approximately 10 bins at 8$\mu$m resolution, 5 bins at 16$\mu$m resolution
#' @param min_cluster_area_um2 Minimum cluster area in $\mu$m$^2$ (default: 1280)
#'   Approximately 20 bins at 8$\mu$m resolution, 5 bins at 16$\mu$m resolution
#'   Default based on 16$\mu$m standard (5 bins = reasonable minimum cluster)
#' @param batch_var Batch variable (default: "sample_id")
#' @param col_x X coordinate column (default: "array_col") 
#' @param col_y Y coordinate column (default: "array_row")
#' @param name Output column prefix (default: "edge_artifact")
#' @param verbose Print progress (default: TRUE)
#' @param keep_intermediate Keep intermediate columns (default: FALSE)
#'
#' @details
#' IMPORTANT: This function uses array_col/array_row (bin indices), NOT
#' pixel coordinates. This is much more memory efficient.
#' 
#' Buffer width and cluster size are specified in physical units (um, um^2)
#' and automatically converted to bins based on resolution:
#' - 8um resolution: 1 bin = 8*8 um = 64 um^2
#' - 16um resolution: 1 bin = 16*16 um = 256 um^2
#' 
#' Default parameters are designed for 16um resolution:
#' - buffer_width_um = 80 um -> 5 bins at 16um, 10 bins at 8um
#' - min_cluster_area_um2 = 1280 um^2 -> 5 bins at 16um, 20 bins at 8um
#' 
#' @examples
#' library(SpatialExperiment)
#' library(S4Vectors)
#' 
#' # 1. Runnable Example: Mock Data (Try this!)
#' 
#' # Create a mock Visium HD dataset (20x20 grid, representing 320x320 um)
#' n_rows <- 20
#' n_cols <- 20
#' n_bins <- n_rows * n_cols
#' coords <- expand.grid(array_row = 1:n_rows, array_col = 1:n_cols)
#' 
#' # Simulate gene counts: Edge artifact (left 2 cols) has low counts
#' counts <- rep(100, n_bins)
#' is_edge <- coords$array_col <= 2
#' counts[is_edge] <- 10 
#' 
#' # Create SpatialExperiment object
#' spe_hd <- SpatialExperiment(
#'     assays = list(counts = matrix(counts, nrow = 1, ncol = n_bins)),
#'     colData = DataFrame(
#'         sample_id = "mock_hd",
#'         in_tissue = rep(TRUE, n_bins),
#'         sum_gene = counts, 
#'         array_row = coords$array_row,
#'         array_col = coords$array_col
#'     )
#' )
#' 
#' # Run detection for 16um resolution
#' # (Physical buffer 80um / 16um bin size = 5 bins. Artifact is 2 bins wide.)
#' spe_hd <- detectEdgeArtifacts_VisiumHD(
#'     spe_hd, 
#'     resolution = "16um",
#'     qc_metric = "sum_gene"
#' )
#' 
#' # Check results
#' table(spe_hd$edge_artifact_edge)
#' 
#' # 2. Illustrative Examples (Concept only)
#' \dontrun{
#' # Assuming 'spe' is a real SpatialExperiment object
#' 
#' # 8um data with defaults
#' # buffer_width_um = 80 -> 10 bins (80 / 8)
#' # min_cluster_area_um2 = 1280 -> 20 bins
#' spe <- detectEdgeArtifacts_VisiumHD(spe, resolution = "8um")
#' 
#' # 16um data with defaults  
#' # buffer_width_um = 80 -> 5 bins (80 / 16)
#' # min_cluster_area_um2 = 1280 -> 5 bins
#' spe <- detectEdgeArtifacts_VisiumHD(spe, resolution = "16um")
#' 
#' # Custom parameters (physical units)
#' spe <- detectEdgeArtifacts_VisiumHD(
#'   spe, 
#'   resolution = "16um",
#'   buffer_width_um = 100,      # 100 $\mu$m buffer
#'   min_cluster_area_um2 = 2000 # 2000 $\mu$m^2 minimum 
#' )
#' }
#'
#' @importFrom scuttle isOutlier
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom S4Vectors DataFrame
#' 
#' @export
detectEdgeArtifacts_VisiumHD <- function(
    spe,
    resolution,  # REQUIRED: "8um" or "16um"
    qc_metric = "sum_gene",
    samples = "sample_id",
    mad_threshold = 3,
    buffer_width_um = 80,         # Physical distance in um (default: ~10 bins at 8um)
    min_cluster_area_um2 = 1280,  # Physical area in um^2 (default: ~5 bins at 16um)
    batch_var = "sample_id",
    col_x = "array_col",         # Bin indices (NOT pixels!)
    col_y = "array_row",
    name = "edge_artifact",
    verbose = TRUE,
    keep_intermediate = FALSE) {

  
  if (!inherits(spe, "SpatialExperiment")) {
    stop("'spe' must be a SpatialExperiment object")
  }
  
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' required for VisiumHD.\n",
         "Install: install.packages('terra')")
  }
  
  if (missing(resolution)) {
    stop("'resolution' is REQUIRED. Specify '8um' or '16um'")
  }
  
  if (!resolution %in% c("8um", "16um")) {
    stop("'resolution' must be '8um' or '16um'. Got: ", resolution)
  }
  
  # Auto-detect or standardize common QC metric names
  if (qc_metric == "sum_gene" && !"sum_gene" %in% colnames(colData(spe))) {
    # Try common alternatives
    if ("detected" %in% colnames(colData(spe))) {
      qc_metric <- "detected"
      if (verbose) message("Using 'detected' as QC metric (sum_gene not found)")
    } else if ("nFeature" %in% colnames(colData(spe))) {
      qc_metric <- "nFeature"
      if (verbose) message("Using 'nFeature' as QC metric (sum_gene not found)")
    } else {
      stop("QC metric 'sum_gene' not found. Available columns: ",
           paste(colnames(colData(spe)), collapse = ", "), "\n",
           "Please specify qc_metric parameter explicitly.")
    }
  }
  
  # Standardize sum to sum_umi if needed
  if ("sum" %in% colnames(colData(spe)) && !"sum_umi" %in% colnames(colData(spe))) {
    colnames(colData(spe))[colnames(colData(spe)) == "sum"] <- "sum_umi"
    if (verbose) message("Renamed 'sum' to 'sum_umi'")
  }
  
  # Check required columns
  required_cols <- c(qc_metric, samples, "in_tissue", col_x, col_y)
  missing_cols <- setdiff(required_cols, colnames(colData(spe)))
  if (length(missing_cols) > 0) {
    stop("Missing columns: ", paste(missing_cols, collapse = ", "), "\n",
         "Available columns: ", paste(colnames(colData(spe)), collapse = ", "))
  }
  
  # Convert physical units to bins based on resolution
  if (resolution == "8um") {
    bin_size_um <- 8
    bin_area_um2 <- 64    # 8*8
    buffer_width_bins <- round(buffer_width_um / bin_size_um)
    min_cluster_size_bins <- ceiling(min_cluster_area_um2 / bin_area_um2)
  } else {  # 16um
    bin_size_um <- 16
    bin_area_um2 <- 256   # 16*16
    buffer_width_bins <- round(buffer_width_um / bin_size_um)
    min_cluster_size_bins <- ceiling(min_cluster_area_um2 / bin_area_um2)
  }
  
  # Get coordinate ranges from data (data-driven boundaries)
  sampleList <- unique(colData(spe)[[samples]])
  names(sampleList) <- sampleList
  
  x_range <- range(colData(spe)[[col_x]], na.rm = TRUE)
  y_range <- range(colData(spe)[[col_y]], na.rm = TRUE)
  
  coord_min_x <- x_range[1]
  coord_max_x <- x_range[2]
  coord_min_y <- y_range[1]
  coord_max_y <- y_range[2]
  
  if (verbose) {
    message("================================================================")
    message("VisiumHD Edge Artifact Detection")
    message("================================================================")
    message(sprintf("Resolution: %s (bin size = %d*%d um)", 
                resolution, bin_size_um, bin_size_um))
    message(sprintf("Coordinate Range: X[%d-%d], Y[%d-%d] (bin indices)", 
                    coord_min_x, coord_max_x, coord_min_y, coord_max_y))
    message(sprintf("Buffer Width: %d um -> %d bins", 
                    buffer_width_um, buffer_width_bins))
    message(sprintf("Min Cluster Area: %d um2 -> %d bins (%.1f bins at 16um)", 
                    min_cluster_area_um2, min_cluster_area_bins, bin_area_um2/256))
    message(sprintf("Samples: %d", length(sampleList)))
  }
  
  if (verbose) message("\n--- STEP 1: Outlier Detection ---")
  
  lg10_metric <- paste0("lg10_", qc_metric)
  colData(spe)[[lg10_metric]] <- log10(colData(spe)[[qc_metric]] + 1)
  
  outlier_binary_col <- paste0(name, "_outlier_binary")
  
  outlier_result <- scuttle::isOutlier(
    colData(spe)[[lg10_metric]],
    subset = colData(spe)$in_tissue,
    batch = colData(spe)[[samples]],
    type = "lower",
    nmads = mad_threshold
  )
  
  thresholds <- attr(outlier_result, "thresholds")
  colData(spe)[[outlier_binary_col]] <- FALSE
  
  for (s in sampleList) {
    sample_cutoff <- thresholds["lower", s]
    s_mask <- colData(spe)[[samples]] == s
    colData(spe)[[outlier_binary_col]][s_mask] <-
      colData(spe)[[lg10_metric]][s_mask] < sample_cutoff
  }
  
  total_outliers <- sum(colData(spe)[[outlier_binary_col]], na.rm = TRUE)
  if (verbose) {
    message(sprintf("  Total outliers detected: %d bins", total_outliers))
  }
  
  if (verbose) {
    message("\n--- STEP 2: Edge Artifact Detection (Buffer Zone) ---")
  }
  
  colData(spe)[[paste0(name, "_edge")]] <- FALSE
  
  for (sample_name in sampleList) {
    sample_mask <- colData(spe)[[samples]] == sample_name
    
    # Get bin indices for this sample
    current_x <- colData(spe)[sample_mask, col_x]
    current_y <- colData(spe)[sample_mask, col_y]
    
    # Define buffer zones near tissue boundaries
    is_left_edge <- current_x <= (coord_min_x + buffer_width_bins)
    is_right_edge <- current_x >= (coord_max_x - buffer_width_bins)
    is_top_edge <- current_y <= (coord_min_y + buffer_width_bins)
    is_bottom_edge <- current_y >= (coord_max_y - buffer_width_bins)
    
    in_buffer_zone <- is_left_edge | is_right_edge | is_top_edge | is_bottom_edge
    
    # Edge artifacts = outliers in buffer zone
    is_edge_artifact <- in_buffer_zone & 
                        colData(spe)[sample_mask, outlier_binary_col]
    
    # Assign results
    global_indices <- which(sample_mask)
    colData(spe)[[paste0(name, "_edge")]][global_indices] <- is_edge_artifact
    
    if (verbose) {
      n_edge <- sum(is_edge_artifact, na.rm = TRUE)
      message(sprintf("  %s: %d edge artifact bins", sample_name, n_edge))
    }
  }
  
  total_edge <- sum(colData(spe)[[paste0(name, "_edge")]], na.rm = TRUE)

  if (verbose) {
    message("\n--- STEP 3: Interior Problem Area Detection ---")
    message("  Method: Morphological processing + clustering")
  }
  
  colData(spe)[[paste0(name, "_problem_id")]] <- NA_character_
  colData(spe)[[paste0(name, "_problem_size")]] <- 0
  
  problem_areas_list <- lapply(sampleList, function(sample_name) {
    
    tmp_dframe <- colData(spe)[colData(spe)[[samples]] == sample_name, ]
    
    # Prepare data: use bin indices directly (NO offset!)
    xyz_df <- data.frame(
      x = tmp_dframe[[col_x]],  # Bin indices
      y = tmp_dframe[[col_y]],
      outlier = tmp_dframe[[outlier_binary_col]],
      stringsAsFactors = FALSE
    )
    rownames(xyz_df) <- rownames(tmp_dframe)
    
    # Morphological clustering
    all_clusters <- problemAreas_WithMorphology_terra(
      .xyz = xyz_df,
      uniqueIdentifier = sample_name,
      min_cluster_size = min_cluster_size_bins,
      resolution = resolution
    )
    
    if (nrow(all_clusters) == 0) {
      if (verbose) {
        message(sprintf("  %s: No problem areas detected", sample_name))
      }
      return(data.frame(
        spotcode = character(0),
        clumpID = character(0),
        clumpSize = numeric(0)
      ))
    }
    
    # Separate edge vs interior clusters based on cluster center
    cluster_ids <- unique(all_clusters$clumpID)
    interior_clusters <- character(0)
    edge_clusters <- character(0)
    
    for (cid in cluster_ids) {
      cluster_spots <- all_clusters$spotcode[all_clusters$clumpID == cid]
      cluster_coords <- tmp_dframe[cluster_spots, c(col_x, col_y), drop = FALSE]
      
      # Cluster center (median is robust)
      center_x <- median(cluster_coords[[col_x]])
      center_y <- median(cluster_coords[[col_y]])
      
      # Check if center is in buffer zone
      in_buffer <- (center_x <= (coord_min_x + buffer_width_bins)) |
                   (center_x >= (coord_max_x - buffer_width_bins)) |
                   (center_y <= (coord_min_y + buffer_width_bins)) |
                   (center_y >= (coord_max_y - buffer_width_bins))
      
      if (in_buffer) {
        edge_clusters <- c(edge_clusters, cid)
      } else {
        interior_clusters <- c(interior_clusters, cid)
      }
    }
    
    # Keep only interior clusters
    interior_problem_areas <- all_clusters[
      all_clusters$clumpID %in% interior_clusters,
      , drop = FALSE
    ]
    
    if (verbose) {
      n_interior_bins <- nrow(interior_problem_areas)
      n_interior_clusters <- length(interior_clusters)
      n_edge_clusters <- length(edge_clusters)
      
      message(sprintf("  %s: %d interior bins in %d clusters (%d edge filtered)",
                      sample_name, n_interior_bins, n_interior_clusters,
                      n_edge_clusters))
    }
    
    return(interior_problem_areas)
  })
  
  # Combine all problem areas
  all_problem_areas <- do.call(rbind, problem_areas_list)
  
  if (nrow(all_problem_areas) > 0) {
    colData(spe)[all_problem_areas$spotcode, paste0(name, "_problem_id")] <-
      all_problem_areas$clumpID
    colData(spe)[all_problem_areas$spotcode, paste0(name, "_problem_size")] <-
      all_problem_areas$clumpSize
  }
  
  total_interior <- sum(!is.na(colData(spe)[[paste0(name, "_problem_id")]]))
  
  if (!keep_intermediate) {
    intermediate_cols <- c(lg10_metric, outlier_binary_col)
    existing_cols <- intersect(intermediate_cols, colnames(colData(spe)))
    
    if (length(existing_cols) > 0) {
      colData(spe) <- colData(spe)[, !colnames(colData(spe)) %in% existing_cols,
                                    drop = FALSE]
    }
  }
  
  
  if (verbose) {
    message("\n================================================================")
    message("Detection Complete!")
    message("================================================================")
    message(sprintf("Total outliers detected:      %d bins", total_outliers))
    message(sprintf("Edge artifacts (in buffer):   %d bins", total_edge))
    message(sprintf("Interior problem areas:       %d bins", total_interior))
    
    if (total_outliers > 0) {
      classification_rate <- 100 * (total_edge + total_interior) / total_outliers
      message(sprintf("Classification rate:          %.1f%%", 
                      classification_rate))
    }
    message("================================================================")
  }
  
  return(spe)
}