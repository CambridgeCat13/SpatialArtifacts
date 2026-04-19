# Unified Edge Artifact Detection

A convenient wrapper that routes to platform-specific edge detection
methods.

## Usage

``` r
detectEdgeArtifacts(spe, platform = c("visium", "visiumhd"), ...)
```

## Arguments

- spe:

  A SpatialExperiment object.

- platform:

  Character string: "visium" or "visiumhd" (case insensitive).

- ...:

  Additional arguments passed to the specific function:

  - For **Visium**: `edge_threshold`, `min_cluster_size`, etc.

  - For **VisiumHD**: `resolution` (REQUIRED), `buffer_width_um`, etc.

## Value

A
[`SpatialExperiment`](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment.html)
object with artifact detection columns (`_edge`, `_problem_id`,
`_problem_size`) added to `colData`.

## Examples

``` r
# Load example data
data(spe_vignette)

# 1. Standard Visium example (runnable)
spe_visium <- detectEdgeArtifacts(spe_vignette, 
platform = "visium", qc_metric = "sum", batch_var = "sample_id")
#> Detecting edges...
#>   Sample V11L05-335_C1: 74 edge spots detected
#> Finding problem areas...
#> Removed intermediate columns: lg10_sum, sum_3MAD_outlier_sample, sum_3MAD_outlier_binary
#> Edge artifact detection completed!
#>   Total edge spots: 74
#>   Total problem area spots: 78
# 2. Visium HD example (wrapped in  to avoid execution without HD data)
# \donttest{
# spe_hd <- detectEdgeArtifacts(spe_hd_example, platform = "visiumhd", resolution = "16um")
# }
```
