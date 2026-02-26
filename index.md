# SpatialArtifacts

The goal of `SpatialArtifacts` is to detect interior and edge artifacts,
such as dry spots caused by incomplete reagent coverage or tissue
handling, in spatial transcriptomics data. The package currently
supports the 10x Genomics `Visium` and `VisiumHD` platforms.

If you experience any issues using the package or would like to make a
suggestion, please open an issue on the [GitHub
repository](https://github.com/CambridgeCat13/SpatialArtifacts/issues).

To find more information, please visit the [documentation
website](https://cambridgecat13.github.io/SpatialArtifacts).

### Key Features

- **Multi-platform support**: Works on both standard **10x Visium** and
  high-resolution **Visium HD**
- **Morphological detection**: Uses raster-based focal transformations
  (fill, outline, star-pattern) to intelligently identify artifact
  clusters
- **Hierarchical classification**: Categorizes artifacts into actionable
  groups (e.g., *Large Edge Artifact*, *Small Interior Artifact*)
- **Fast and efficient**: Optimized with the
  [`terra`](https://CRAN.R-project.org/package=terra) R package for
  handling large datasets

## Installation

You can install the latest version of `SpatialArtifacts` from
Bioconductor with the following code:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("SpatialArtifacts")
```

You can install the development version of SpatialArtifacts from
[GitHub](https://github.com/CambridgeCat13/SpatialArtifacts) with:

``` r
# install.packages("pak")
pak::pak("CambridgeCat13/SpatialArtifacts")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(SpatialArtifacts)
library(SpatialExperiment)

# 1. Detect artifacts
# Option A: Standard Visium (Hexagonal Grid)
# Just specify platform = "visium"
spe <- detectEdgeArtifacts(spe, platform = "visium", qc_metric = "sum_umi")

# Option B: Visium HD (Square Grid)
# Specify platform = "visiumhd" AND the required resolution ("8um" or "16um")
# spe <- detectEdgeArtifacts(spe, platform = "visiumhd", resolution = "8um")

# 2. Classify results (Platform independent)
# Note: For Visium HD, remember to increase min_spots (e.g., min_spots = 400)
spe <- classifyEdgeArtifacts(spe)

# 3. View classification
table(spe$edge_artifact_classification)
```

## Tutorials

A detailed tutorial is available in the package vignette from
Bioconductor. A direct link to the tutorial / package vignette is
available
[here](https://cambridgecat13.github.io/SpatialArtifacts/articles/hippocampus-edge-detection.html).

## Development tools

- Continuous code testing is possible thanks to GitHub actions.
- The [documentation
  website](https://cambridgecat13.github.io/SpatialArtifacts/) is
  automatically updated thanks to `BiocStyle::CRANpkg('pkgdown')`.
- The documentation is formatted thanks to
  `BiocStyle::CRANpkg('devtools')` and `BiocStyle::CRANpkg('roxygen2')`.
- This package was developed using `BiocStyle::Biocpkg('biocthis')`.

## Contributors

- Harriet Jiali He
- Jacqueline Thompson
- Michael Totty
- [Stephanie C. Hicks](https://stephaniehicks.com)
