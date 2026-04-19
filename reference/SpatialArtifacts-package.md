# SpatialArtifacts: Detect and Classify Spatial Artifacts

`SpatialArtifacts` provides a robust, data-driven two-step workflow to
identify, classify, and handle spatial artifacts in spatial
transcriptomics data across multiple platforms including 10x Genomics
Visium (Standard and HD). It combines median absolute deviation
(MAD)-based outlier detection with morphological image processing to
flag problematic regions such as edge dryspots and interior artifacts
caused by incomplete reagent coverage.

## Main Functions

- [`detectEdgeArtifacts`](https://cambridgecat13.github.io/SpatialArtifacts/reference/detectEdgeArtifacts.md):
  The primary wrapper function to detect potential artifact spots.
  Automatically routes to platform-specific methods based on the
  `platform` argument. Outputs three columns to `colData`: `*_edge`,
  `*_problem_id`, and `*_problem_size`.

- [`classifyEdgeArtifacts`](https://cambridgecat13.github.io/SpatialArtifacts/reference/classifyEdgeArtifacts.md):
  Hierarchically classifies detected artifacts by location (edge vs.
  interior) and size (large vs. small). Outputs a single
  `*_classification` column.

## Typical Workflow

    # Step 1: Detect artifacts
    spe <- detectEdgeArtifacts(spe, platform = "visium", qc_metric = "sum_gene")

    # Step 2: Classify artifacts
    spe <- classifyEdgeArtifacts(spe, min_spots = 20)

## Platform-Specific Usage

[`detectEdgeArtifacts`](https://cambridgecat13.github.io/SpatialArtifacts/reference/detectEdgeArtifacts.md)
requires users to specify their platform:

- **Standard Visium** (`platform = "visium"`): Uses hexagonal grid
  layout. The default `shifted = FALSE` is correct for standard Space
  Ranger `array_row`/`array_col` outputs.

- **Visium HD** (`platform = "visiumhd"`): Uses square grid layout.
  Requires the `resolution` parameter (`"16um"` or `"8um"`). Parameters
  are specified in physical units (micrometers).

## Input Data

All functions accept a
[`SpatialExperiment`](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment.html)
object. QC metrics (e.g., library size, detected genes) should be
precomputed, for example using
[`addPerCellQCMetrics`](https://rdrr.io/pkg/scuttle/man/addPerCellQCMetrics.html).

## See also

Useful links:

- <https://github.com/CambridgeCat13/SpatialArtifacts>

- Report bugs at
  <https://github.com/CambridgeCat13/SpatialArtifacts/issues>

## Author

**Maintainer**: Harriet Jiali He <jhe46@jh.edu>
([ORCID](https://orcid.org/0009-0003-7827-2735))

Authors:

- Jacqueline R. Thompson <jthom338@jh.edu>

- Michael Totty <mtotty2@jh.edu>

- Stephanie C. Hicks <shicks19@jhu.edu>
  ([ORCID](https://orcid.org/0000-0002-7858-0231)) \[funder\]
