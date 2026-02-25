#' Unified Edge Artifact Detection
#'
#' A convenient wrapper that routes to platform-specific edge detection methods.
#' 
#' @param spe A SpatialExperiment object
#' @param platform Character string: "visium" or "visiumhd" (case insensitive)
#' @param ... Additional arguments passed to the specific function:
#'   \itemize{
#'     \item For **Visium**: \code{edge_threshold}, \code{min_cluster_size}, etc.
#'     \item For **VisiumHD**: \code{resolution} (REQUIRED), \code{buffer_width_um}, etc.
#'   }
#'
#' @return A SpatialExperiment object with artifact annotations.
#'
#' @examples
#' # 1. Standard Visium
#' # spe <- detectEdgeArtifacts(spe, platform = "Visium")
#'
#' # 2. Visium HD (Must provide resolution)
#' # spe <- detectEdgeArtifacts(spe, platform = "VisiumHD", resolution = "16um")
#'
#' @export
detectEdgeArtifacts <- function(spe, platform = c("visium", "visiumhd"), ...) {
  platform <- match.arg(tolower(platform), c("visium", "visiumhd"))
  
  if (platform == "visium") {
    return(detectEdgeArtifacts_Visium(spe, ...))
    
  } else if (platform == "visiumhd") {
    args <- list(...)
    if (is.null(args$resolution)) {
      stop("Error: 'resolution' is REQUIRED when platform = 'visiumhd'.\n",
           "Usage: detectEdgeArtifacts(spe, platform='visiumhd', resolution='16um')")
    }

    return(detectEdgeArtifacts_VisiumHD(spe, ...))
  }
}