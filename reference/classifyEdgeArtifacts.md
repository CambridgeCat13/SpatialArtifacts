# Classify Edge Artifacts into Hierarchical Categories

Classifies detected artifacts based on their location (edge vs.
interior) and size (large vs. small). This function works downstream of
[`detectEdgeArtifacts()`](https://cambridgecat13.github.io/SpatialArtifacts/reference/detectEdgeArtifacts.md).

## Usage

``` r
classifyEdgeArtifacts(
  spe,
  qc_metric = "sum_umi",
  samples = "sample_id",
  min_spots = 20,
  name = "edge_artifact",
  exclude_slides = NULL
)
```

## Arguments

- spe:

  A SpatialExperiment object that has been processed with
  [`detectEdgeArtifacts()`](https://cambridgecat13.github.io/SpatialArtifacts/reference/detectEdgeArtifacts.md).

- qc_metric:

  Character string specifying the QC metric column used for validation
  (default: "sum_umi"). Note: This column must exist, but is not
  directly used for classification logic.

- samples:

  Character string specifying the sample ID column name (default:
  "sample_id").

- min_spots:

  Minimum number of spots for an artifact to be classified as "large"
  (default: 20).

- name:

  Character string matching the `name` argument used in
  [`detectEdgeArtifacts()`](https://cambridgecat13.github.io/SpatialArtifacts/reference/detectEdgeArtifacts.md)
  (default: "edge_artifact").

- exclude_slides:

  Character vector of slide IDs to exclude from edge detection (default:
  NULL). Spots on these slides will be forced to FALSE for edge
  artifacts.

## Value

A SpatialExperiment object with additional classification columns in
`colData`:

- `[name]_true_edges`: Logical, indicating edge artifacts after applying
  slide exclusions.

- `[name]_classification`: Character, containing hierarchical
  categories:

  - `"not_artifact"`: Normal, high-quality spots.

  - `"large_edge_artifact"`: Clusters \> min_spots touching the slide
    boundary.

  - `"small_edge_artifact"`: Clusters \<= min_spots touching the slide
    boundary.

  - `"large_interior_artifact"`: Clusters \> min_spots in the tissue
    interior.

  - `"small_interior_artifact"`: Clusters \<= min_spots in the tissue
    interior.

## Details

**Parameter Recommendations:**

The `min_spots` threshold should scale with platform resolution to
represent similar physical artifact sizes:

- **Standard Visium (55µm bins):** `min_spots = 20-40`

  - Physical area: ~0.06-0.12 mm²

  - Rationale: Artifacts \<20 spots likely represent isolated noise

  - Typical edge artifacts: 50-200 spots

- **VisiumHD 16µm bins:** `min_spots = 100-200`

  - Physical area: ~0.026-0.051 mm² (comparable to Visium)

  - Rationale: Higher density requires proportionally higher threshold

  - Typical edge artifacts: 500-2000 bins

- **VisiumHD 8µm bins:** `min_spots = 400-800`

  - Physical area: ~0.026-0.051 mm² (comparable to Visium)

  - Rationale: 4× density of 16µm bins requires 4× threshold

  - Typical edge artifacts: 2000-8000 bins

**Practical Guideline:** For VisiumHD data, start with:
`min_spots = 20 × (55µm / bin_size)²`

This formula maintains constant physical area thresholds across
resolutions.

**Context-Specific Tuning:**

- Increase `min_spots` for noisy data (prevents over-flagging small
  clusters)

- Decrease `min_spots` for high-quality data (captures smaller genuine
  artifacts)

- Visualize intermediate results to validate threshold appropriateness

## See also

[`detectEdgeArtifacts`](https://cambridgecat13.github.io/SpatialArtifacts/reference/detectEdgeArtifacts.md)

## Examples

``` r
library(SpatialExperiment)
#> Loading required package: SingleCellExperiment
#> Loading required package: SummarizedExperiment
#> Loading required package: MatrixGenerics
#> Loading required package: matrixStats
#> 
#> Attaching package: ‘MatrixGenerics’
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
#>     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
#>     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
#>     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
#>     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
#>     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
#>     colWeightedMeans, colWeightedMedians, colWeightedSds,
#>     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
#>     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
#>     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
#>     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
#>     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
#>     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
#>     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
#>     rowWeightedSds, rowWeightedVars
#> Loading required package: GenomicRanges
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Loading required package: generics
#> 
#> Attaching package: ‘generics’
#> The following objects are masked from ‘package:base’:
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: ‘BiocGenerics’
#> The following objects are masked from ‘package:stats’:
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from ‘package:base’:
#> 
#>     Filter, Find, Map, Position, Reduce, anyDuplicated, aperm, append,
#>     as.data.frame, basename, cbind, colnames, dirname, do.call,
#>     duplicated, eval, evalq, get, grep, grepl, is.unsorted, lapply,
#>     mapply, match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
#>     rank, rbind, rownames, sapply, saveRDS, table, tapply, unique,
#>     unsplit, which.max, which.min
#> Loading required package: S4Vectors
#> 
#> Attaching package: ‘S4Vectors’
#> The following object is masked from ‘package:utils’:
#> 
#>     findMatches
#> The following objects are masked from ‘package:base’:
#> 
#>     I, expand.grid, unname
#> Loading required package: IRanges
#> Loading required package: Seqinfo
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 
#> Attaching package: ‘Biobase’
#> The following object is masked from ‘package:MatrixGenerics’:
#> 
#>     rowMedians
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     anyMissing, rowMedians
library(S4Vectors)

# --- Create a Mock SPE with "Detected" Results ---
# We simulate the output that detectEdgeArtifacts() would produce
n_spots <- 5
spe <- SpatialExperiment(
    colData = DataFrame(
        sample_id = "sample01",
        sum_umi = rep(100, n_spots), 
        edge_artifact_edge = c(TRUE, TRUE, FALSE, FALSE, FALSE),
        edge_artifact_problem_id = c("Cluster1", "Cluster1", "Cluster2", "Cluster3", NA),
        edge_artifact_problem_size = c(50, 50, 10, 30, 0)
    )
)

# Review the mock scenarios:
# Spot 1: Edge=TRUE, Size=50 -> Expect "large_edge_artifact"
# Spot 2: Edge=TRUE, Size=50 -> Expect "large_edge_artifact"
# Spot 3: Edge=FALSE, Size=10 -> Expect "small_interior_artifact"
# Spot 4: Edge=FALSE, Size=30 -> Expect "large_interior_artifact"
# Spot 5: Edge=FALSE, Size=0  -> Expect "not_artifact"

# --- Run Classification ---
# Set threshold to 20 to separate large/small
spe <- classifyEdgeArtifacts(spe, min_spots = 20, name = "edge_artifact")
#> Classifying artifacts spots...
#> Classification added: edge_artifact_classification
#> 
#> Classification summary:
#>   large_edge_artifact: 2 spots
#>   large_interior_artifact: 1 spots
#>   not_artifact: 1 spots
#>   small_interior_artifact: 1 spots
table(spe$edge_artifact_classification)
#> 
#>     large_edge_artifact large_interior_artifact            not_artifact 
#>                       2                       1                       1 
#> small_interior_artifact 
#>                       1 
colData(spe)[, c("edge_artifact_problem_size", "edge_artifact_classification")]
#> DataFrame with 5 rows and 2 columns
#>   edge_artifact_problem_size edge_artifact_classification
#>                    <numeric>                  <character>
#> 1                         50          large_edge_artifact
#> 2                         50          large_edge_artifact
#> 3                         10       small_interior_artif..
#> 4                         30       large_interior_artif..
#> 5                          0                 not_artifact
```
