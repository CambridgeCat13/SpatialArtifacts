# Detect edge artifacts in spatial transcriptomics data

This function identifies edge artifacts and problem areas in spatial
transcriptomics data by analyzing QC metrics and spatial patterns.

## Usage

``` r
detectEdgeArtifacts_Visium(
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
  keep_intermediate = FALSE
)
```

## Arguments

- spe:

  A SpatialExperiment object containing spatial transcriptomics data

- qc_metric:

  Character string specifying the QC metric column name to analyze
  (default: "sum_gene")

- samples:

  Character string specifying the sample ID column name (default:
  "sample_id")

- mad_threshold:

  Numeric value for MAD threshold for outlier detection (default: 3)

- edge_threshold:

  Numeric threshold for edge detection (default: 0.75)

- min_cluster_size:

  Minimum cluster size for morphological cleaning (default: 40)

- shifted:

  Logical indicating whether to apply coordinate adjustment for
  hexagonal arrays (default: FALSE)

- batch_var:

  Character specifying batch variable for outlier detection ("slide",
  "sample_id", or "both", default: "both")

- name:

  Character string for naming output columns (default: "edge_artifact")

- verbose:

  Logical indicating whether to print progress messages (default: TRUE)

- keep_intermediate:

  Logical indicating whether to keep intermediate outlier detection
  columns (default: FALSE)

## Value

A SpatialExperiment object with additional columns in colData:

- [terra::name](https://rspatial.github.io/terra/reference/names.html)\_edge:

  Logical indicating spots identified as edges

- [terra::name](https://rspatial.github.io/terra/reference/names.html)\_problem_id:

  Character identifying problem area clusters

- [terra::name](https://rspatial.github.io/terra/reference/names.html)\_problem_size:

  Numeric size of problem area clusters

## Examples

``` r
library(SummarizedExperiment)
library(SpatialExperiment)
library(S4Vectors)

# Create a minimal mock SpatialExperiment (4x4 grid with edge artifact)
set.seed(123)
n_spots <- 16
coords <- expand.grid(row = 1:4, col = 1:4)

# Simulate counts with lower values at edges (top row)
mock_counts <- rpois(n_spots, lambda = 500)
mock_counts[1:4] <- rpois(4, lambda = 50) # Edge artifact

spe_mock <- SpatialExperiment::SpatialExperiment(
  assays = list(counts = matrix(rpois(n_spots * 10, lambda = 5),
    nrow = 10, ncol = n_spots
  )),
  colData = DataFrame(
    in_tissue = rep(TRUE, n_spots),
    sum_gene = mock_counts,
    sum_umi = mock_counts, # Add sum_umi for classify function
    sample_id = "mock_sample",
    slide = "mock_slide",
    array_row = coords$row,
    array_col = coords$col
  ),
  spatialCoords = as.matrix(coords)
)

colnames(spe_mock) <- paste0("spot_", seq_len(n_spots))
rownames(spe_mock) <- paste0("gene_", 1:10)

# Detect edge artifacts
spe_detected <- detectEdgeArtifacts_Visium(
  spe_mock,
  qc_metric = "sum_gene",
  samples = "sample_id",
  mad_threshold = 3,
  min_cluster_size = 1,
  name = "edge_artifact"
)
#> Detecting edges...
#>   Sample mock_sample: 4 edge spots detected
#> Finding problem areas...
#> Removed intermediate columns: lg10_sum_gene, sum_gene_3MAD_outlier_slide, sum_gene_3MAD_outlier_sample, sum_gene_3MAD_outlier_binary
#> Edge artifact detection completed!
#>   Total edge spots: 4
#>   Total problem area spots: 4

# Check detection results
table(spe_detected$edge_artifact_edge)
#> 
#> FALSE  TRUE 
#>    12     4 
head(colData(spe_detected)[, c(
  "edge_artifact_edge",
  "edge_artifact_problem_id"
)])
#> DataFrame with 6 rows and 2 columns
#>        edge_artifact_edge edge_artifact_problem_id
#>                 <logical>              <character>
#> spot_1               TRUE            mock_sample_1
#> spot_2               TRUE            mock_sample_1
#> spot_3               TRUE            mock_sample_1
#> spot_4               TRUE            mock_sample_1
#> spot_5              FALSE                       NA
#> spot_6              FALSE                       NA
```
