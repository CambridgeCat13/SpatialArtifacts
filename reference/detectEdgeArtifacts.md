# Unified Edge Artifact Detection

A convenient wrapper that routes to platform-specific edge detection
methods.

## Usage

``` r
detectEdgeArtifacts(spe, platform = c("visium", "visiumhd"), ...)
```

## Arguments

- spe:

  A SpatialExperiment object

- platform:

  Character string: "visium" or "visiumhd" (case insensitive)

- ...:

  Additional arguments passed to the specific function:

  - For **Visium**: `edge_threshold`, `min_cluster_size`, etc.

  - For **VisiumHD**: `resolution` (REQUIRED), `buffer_width_um`, etc.

## Value

A SpatialExperiment object with artifact annotations.

## Examples

``` r
# 1. Standard Visium
# spe <- detectEdgeArtifacts(spe, platform = "Visium")

# 2. Visium HD (Must provide resolution)
# spe <- detectEdgeArtifacts(spe, platform = "VisiumHD", resolution = "16um")
```
