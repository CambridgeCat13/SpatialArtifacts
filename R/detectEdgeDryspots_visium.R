#' Detect edge artifacts in spatial transcriptomics data
#'
#' This function identifies edge artifacts and problem areas in spatial
#' transcriptomics data by analyzing QC metrics and spatial patterns.
#'
#' @param spe A SpatialExperiment object containing spatial transcriptomics data
#' @param qc_metric Character string specifying the QC metric column name to analyze (default: "sum_gene")
#' @param samples Character string specifying the sample ID column name (default: "sample_id")
#' @param mad_threshold Numeric value for MAD threshold for outlier detection (default: 3)
#' @param edge_threshold Numeric threshold for edge detection (default: 0.75)
#' @param min_cluster_size Minimum cluster size for morphological cleaning (default: 40)
#' @param shifted Logical indicating whether to apply coordinate adjustment for hexagonal arrays (default: FALSE)
#' @param batch_var Character specifying batch variable for outlier detection ("slide", "sample_id", or "both", default: "both")
#' @param name Character string for naming output columns (default: "edge_artifact")
#' @param verbose Logical indicating whether to print progress messages (default: TRUE)
#' @param keep_intermediate Logical indicating whether to keep intermediate outlier detection columns (default: FALSE)
#'
#' @return A SpatialExperiment object with additional columns in colData:
#'   \item{[name]_edge}{Logical indicating spots identified as edges}
#'   \item{[name]_problem_id}{Character identifying problem area clusters}
#'   \item{[name]_problem_size}{Numeric size of problem area clusters}
#'
#' @examples
#' library(SummarizedExperiment)
#' library(SpatialExperiment)
#' library(S4Vectors)
#'
#' # Create a minimal mock SpatialExperiment (4x4 grid with edge artifact)
#' set.seed(123)
#' n_spots <- 16
#' coords <- expand.grid(row = 1:4, col = 1:4)
#'
#' # Simulate counts with lower values at edges (top row)
#' mock_counts <- rpois(n_spots, lambda = 500)
#' mock_counts[1:4] <- rpois(4, lambda = 50)  # Edge artifact
#'
#' spe_mock <- SpatialExperiment::SpatialExperiment(
#'   assays = list(counts = matrix(rpois(n_spots * 10, lambda = 5),
#'                                  nrow = 10, ncol = n_spots)),
#'   colData = DataFrame(
#'     in_tissue = rep(TRUE, n_spots),
#'     sum = mock_counts,
#'     sum_umi = mock_counts, # Add sum_umi for classify function
#'     sample_id = "mock_sample",
#'     slide = "mock_slide",
#'     array_row = coords$row,
#'     array_col = coords$col
#'   ),
#'   spatialCoords = as.matrix(coords)
#' )
#'
#' colnames(spe_mock) <- paste0("spot_", 1:n_spots)
#' rownames(spe_mock) <- paste0("gene_", 1:10)
#'
#' # Detect edge artifacts
#' spe_detected <- detectEdgeArtifacts_Visium(
#'   spe_mock,
#'   qc_metric = "sum_gene",
#'   samples = "sample_id",
#'   mad_threshold = 3,
#'   min_cluster_size = 1,
#'   name = "edge_artifact"
#' )
#'
#' # Check detection results
#' table(spe_detected$edge_artifact_edge)
#' head(colData(spe_detected)[, c("edge_artifact_edge",
#'                                 "edge_artifact_problem_id")])
#' @import SpatialExperiment
#' @importFrom dplyr %>% sym
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom scuttle isOutlier
#' 
#' @export
detectEdgeArtifacts_Visium <- function(
    spe,
    qc_metric = "sum_gene",
    samples = "sample_id",
    mad_threshold = 3,
    edge_threshold = 0.75,
    min_cluster_size = 40,
    shifted = FALSE,
    batch_var = "both",
    name = "edge_artifact", 
    verbose = TRUE,
    keep_intermediate = FALSE) {

  lg10_metric <- paste0("lg10_", qc_metric)
  colData(spe)[[lg10_metric]] <- log10(colData(spe)[[qc_metric]] + 1)

  if (batch_var %in% c("slide", "both")) {
    if ("slide" %in% colnames(colData(spe))) {
      outlier_slide_col <- paste0(qc_metric, "_", mad_threshold, "MAD_outlier_slide")
      colData(spe)[[outlier_slide_col]] <- isOutlier(colData(spe)[[lg10_metric]], subset = colData(spe)$in_tissue, batch = colData(spe)$slide, type = "lower", nmads = mad_threshold)
    }
  }
  if (batch_var %in% c("sample_id", "both")) {
    outlier_sample_col <- paste0(qc_metric, "_", mad_threshold, "MAD_outlier_sample")
    colData(spe)[[outlier_sample_col]] <- isOutlier(colData(spe)[[lg10_metric]], subset = colData(spe)$in_tissue, batch = colData(spe)[[samples]], type = "lower", nmads = mad_threshold)
  }

  outlier_binary_col <- paste0(qc_metric, "_", mad_threshold, "MAD_outlier_binary")
  if (batch_var == "both") {
    colData(spe)[[outlier_binary_col]] <- colData(spe)[[paste0(qc_metric, "_", mad_threshold, "MAD_outlier_slide")]] | colData(spe)[[paste0(qc_metric, "_", mad_threshold, "MAD_outlier_sample")]]
  } else if (batch_var == "slide") {
    colData(spe)[[outlier_binary_col]] <- colData(spe)[[paste0(qc_metric, "_", mad_threshold, "MAD_outlier_slide")]]
  } else {
    colData(spe)[[outlier_binary_col]] <- colData(spe)[[paste0(qc_metric, "_", mad_threshold, "MAD_outlier_sample")]]
  }
  colData(spe)[[outlier_binary_col]][!colData(spe)$in_tissue] <- FALSE

  sampleList <- unique(colData(spe)[[samples]])
  names(sampleList) <- sampleList

  if (verbose) message("Detecting edges...")

  # Efficiently calculate edges for all samples at once
  edge_spots_list <- lapply(sampleList, function(sample_name) {
    tmp_dframe <- colData(spe)[colData(spe)[[samples]] == sample_name, ]
    xyz_df <- as.data.frame(tmp_dframe[, c("array_row", "array_col", outlier_binary_col)])
    rownames(xyz_df) <- rownames(tmp_dframe)
    
    clumpEdges(
      xyz_df,
      offTissue = rownames(tmp_dframe)[tmp_dframe$in_tissue == FALSE],
      shifted = shifted,
      edge_threshold = edge_threshold,
      min_cluster_size = min_cluster_size
    )
  })
  
  # Initialize the result column
  colData(spe)[[paste0(name, "_edge")]] <- FALSE
  # Accurately assign results back to the spe object
  for(sample_name in names(edge_spots_list)) {
    spots_to_flag <- edge_spots_list[[sample_name]]
    if (length(spots_to_flag) > 0) {
      colData(spe)[spots_to_flag, paste0(name, "_edge")] <- TRUE
      if (verbose) message(sprintf("  Sample %s: %d edge spots detected", sample_name, length(spots_to_flag)))
    }
  }

  if (verbose) message("Finding problem areas...")

  # Efficiently find problem areas for all samples
  problem_areas_list <- lapply(sampleList, function(sample_name) {
    tmp_dframe <- colData(spe)[colData(spe)[[samples]] == sample_name, ]
    xyz_df <- as.data.frame(tmp_dframe[, c("array_row", "array_col", outlier_binary_col)])
    rownames(xyz_df) <- rownames(tmp_dframe)

    problemAreas(
      xyz_df,
      offTissue = rownames(tmp_dframe)[tmp_dframe$in_tissue == FALSE],
      uniqueIdentifier = sample_name,
      shifted = shifted,
      min_cluster_size = min_cluster_size
    )
  })

  # Combine all problem area data frames into one
  all_problem_areas <- do.call(rbind, problem_areas_list)

  # Initialize result columns
  colData(spe)[[paste0(name, "_problem_id")]] <- NA
  colData(spe)[[paste0(name, "_problem_size")]] <- 0
  
  # Accurately assign results in one vectorized operation
  if (nrow(all_problem_areas) > 0) {
    # Match by barcode (spotcode) which should be unique within the problem area dataframe
    colData(spe)[all_problem_areas$spotcode, paste0(name, "_problem_id")] <- all_problem_areas$clumpID
    colData(spe)[all_problem_areas$spotcode, paste0(name, "_problem_size")] <- all_problem_areas$clumpSize
  }
  
  # Clean up intermediate columns if requested
  if (!keep_intermediate) {
    intermediate_cols <- c(
      lg10_metric,
      paste0(qc_metric, "_", mad_threshold, "MAD_outlier_slide"),
      paste0(qc_metric, "_", mad_threshold, "MAD_outlier_sample"),
      outlier_binary_col
    )
    
    # Only remove columns that exist
    existing_cols <- intersect(intermediate_cols, colnames(colData(spe)))
    if (length(existing_cols) > 0) {
      colData(spe) <- colData(spe)[, !colnames(colData(spe)) %in% existing_cols]
      if (verbose) message("Removed intermediate columns: ", paste(existing_cols, collapse = ", "))
    }
  }
  
  if (verbose) {
    total_edge_spots <- sum(colData(spe)[[paste0(name, "_edge")]], na.rm = TRUE)
    total_problem_spots <- sum(!is.na(colData(spe)[[paste0(name, "_problem_id")]]))
    message(sprintf("Edge artifact detection completed!\n  Total edge spots: %d\n  Total problem area spots: %d", 
                    total_edge_spots, total_problem_spots))
  }
  
  return(spe)
}