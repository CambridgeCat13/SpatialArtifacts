#' SpatialArtifacts: Detect and Classify Spatial Artifacts
#'
#' @description
#' `SpatialArtifacts` provides a robust, data-driven two-step workflow to
#' identify, classify, and handle spatial artifacts in spatial transcriptomics
#' data across multiple platforms including 10x Genomics Visium (Standard and
#' HD). It combines median absolute deviation (MAD)-based outlier detection
#' with morphological image processing to flag problematic regions such as
#' edge dryspots and interior artifacts caused by incomplete reagent coverage.
#'
#' @section Main Functions:
#' \itemize{
#'   \item \code{\link{detectEdgeArtifacts}}: The primary wrapper function to
#'     detect potential artifact spots. Automatically routes to
#'     platform-specific methods based on the \code{platform} argument.
#'     Outputs three columns to \code{colData}: \code{*_edge},
#'     \code{*_problem_id}, and \code{*_problem_size}.
#'   \item \code{\link{classifyEdgeArtifacts}}: Hierarchically classifies
#'     detected artifacts by location (edge vs. interior) and size (large
#'     vs. small). Outputs a single \code{*_classification} column.
#' }
#'
#' @section Typical Workflow:
#' \preformatted{
#' # Step 1: Detect artifacts
#' spe <- detectEdgeArtifacts(spe, platform = "visium", qc_metric = "sum_gene")
#'
#' # Step 2: Classify artifacts
#' spe <- classifyEdgeArtifacts(spe, min_spots = 20)
#' }
#'
#' @section Platform-Specific Usage:
#' \code{\link{detectEdgeArtifacts}} requires users to specify their platform:
#' \itemize{
#'   \item \strong{Standard Visium} (\code{platform = "visium"}): Uses
#'     hexagonal grid layout. The default \code{shifted = FALSE} is correct
#'     for standard Space Ranger \code{array_row}/\code{array_col} outputs.
#'   \item \strong{Visium HD} (\code{platform = "visiumhd"}): Uses square
#'     grid layout. Requires the \code{resolution} parameter
#'     (\code{"16um"} or \code{"8um"}). Parameters are specified in physical
#'     units (micrometers).
#' }
#'
#' @section Input Data:
#' All functions accept a \code{\link[SpatialExperiment]{SpatialExperiment}}
#' object. QC metrics (e.g., library size, detected genes) should be
#' precomputed, for example using
#' \code{\link[scuttle]{addPerCellQCMetrics}}.
#'
#' @seealso
#' Useful links:
#' \itemize{
#'   \item \url{https://github.com/CambridgeCat13/SpatialArtifacts}
#'   \item Report bugs at \url{https://github.com/CambridgeCat13/SpatialArtifacts/issues}
#' }
#'
#' @name SpatialArtifacts-package
#' @aliases SpatialArtifacts
"_PACKAGE"