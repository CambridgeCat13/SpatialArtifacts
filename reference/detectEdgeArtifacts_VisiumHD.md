# Detect Edge Artifacts in VisiumHD Data

Detect Edge Artifacts in VisiumHD Data

## Usage

``` r
detectEdgeArtifacts_VisiumHD(
  spe,
  resolution,
  qc_metric = "sum_gene",
  samples = "sample_id",
  mad_threshold = 3,
  buffer_width_um = 80,
  min_cluster_area_um2 = 1280,
  batch_var = "sample_id",
  col_x = "array_col",
  col_y = "array_row",
  name = "edge_artifact",
  verbose = TRUE,
  keep_intermediate = FALSE
)
```

## Arguments

- spe:

  SpatialExperiment object

- resolution:

  Resolution: "8um" or "16um" (REQUIRED)

- qc_metric:

  QC metric column (default: "sum_gene")

- samples:

  Sample ID column (default: "sample_id")

- mad_threshold:

  MAD threshold (default: 3)

- buffer_width_um:

  Buffer zone width in micrometers (default: 80) Approximately 10 bins
  at 8\$u\$m resolution, 5 bins at 16\$u\$m resolution

- min_cluster_area_um2:

  Minimum cluster area in \$u\$m\$^2\$ (default: 1280) Approximately 20
  bins at 8\$u\$m resolution, 5 bins at 16\$u\$m resolution Default
  based on 16\$u\$m standard (5 bins = reasonable minimum cluster)

- batch_var:

  Batch variable (default: "sample_id")

- col_x:

  X coordinate column (default: "array_col")

- col_y:

  Y coordinate column (default: "array_row")

- name:

  Output column prefix (default: "edge_artifact")

- verbose:

  Print progress (default: TRUE)

- keep_intermediate:

  Keep intermediate columns (default: FALSE)

## Details

IMPORTANT: This function uses array_col/array_row (bin indices), NOT
pixel coordinates. This is much more memory efficient.

Buffer width and cluster size are specified in physical units (um, um^2)
and automatically converted to bins based on resolution:

- 8um resolution: 1 bin = 8\*8 um = 64 um^2

- 16um resolution: 1 bin = 16\*16 um = 256 um^2

Default parameters are designed for 16um resolution:

- buffer_width_um = 80 um -\> 5 bins at 16um, 10 bins at 8um

- min_cluster_area_um2 = 1280 um^2 -\> 5 bins at 16um, 20 bins at 8um

## Examples

``` r
library(SpatialExperiment)
library(S4Vectors)

# 1. Runnable Example: Mock Data (Try this!)

# Create a mock Visium HD dataset (20x20 grid, representing 320x320 um)
n_rows <- 20
n_cols <- 20
n_bins <- n_rows * n_cols
coords <- expand.grid(array_row = 1:n_rows, array_col = 1:n_cols)

# Simulate gene counts: Edge artifact (left 2 cols) has low counts
counts <- rep(100, n_bins)
is_edge <- coords$array_col <= 2
counts[is_edge] <- 10 

# Create SpatialExperiment object
spe_hd <- SpatialExperiment(
    assays = list(counts = matrix(counts, nrow = 1, ncol = n_bins)),
    colData = DataFrame(
        sample_id = "mock_hd",
        in_tissue = rep(TRUE, n_bins),
        sum_gene = counts, 
        array_row = coords$array_row,
        array_col = coords$array_col
    )
)

# Run detection for 16um resolution
# (Physical buffer 80um / 16um bin size = 5 bins. Artifact is 2 bins wide.)

colnames(spe_hd) <- paste0("spot_", seq_len(ncol(spe_hd)))
rownames(spe_hd) <- paste0("gene_", seq_len(nrow(spe_hd)))

spe_hd <- detectEdgeArtifacts_VisiumHD(
    spe_hd, 
    resolution = "16um",
    qc_metric = "sum_gene"
)
#> ================================================================
#> VisiumHD Edge Artifact Detection
#> ================================================================
#> Resolution: 16um (bin size = 16*16 um)
#> Coordinate Range: X[1-20], Y[1-20] (bin indices)
#> Buffer Width: 80 um -> 5 bins
#> Min Cluster Area: 1280 um2 -> 5 bins
#> Samples: 1
#> 
#> --- STEP 1: Outlier Detection ---
#>   Total outliers detected: 40 bins
#> 
#> --- STEP 2: Edge Artifact Detection (Buffer Zone) ---
#>   mock_hd: 40 edge artifact bins
#> 
#> --- STEP 3: Interior Problem Area Detection ---
#>   Method: Morphological processing + clustering
#>     Clustering 40 outlier bins...
#>     Creating raster: 20 rows * 20 cols
#>     Applying morphological transformations...
#>     [16um Mode] Adjusted min_cluster_size: 5 -> 2
#>     Clustering outlier regions...
#>     Found 40 clustered cells
#>     Identified 1 unique clusters
#>     Kept 1 clusters (>=2 bins): 40 bins total
#>   mock_hd: 0 interior bins in 0 clusters (1 edge filtered)
#> 
#> ================================================================
#> Detection Complete!
#> ================================================================
#> Total outliers detected:      40 bins
#> Edge artifacts (in buffer):   40 bins
#> Interior problem areas:       0 bins
#> Classification rate:          100.0%
#> ================================================================

# Check results
table(spe_hd$edge_artifact_edge)
#> 
#> FALSE  TRUE 
#>   360    40 

# 2. Illustrative Examples (Concept only)
if (FALSE) { # \dontrun{
# Assuming 'spe' is a real SpatialExperiment object

# 8um data with defaults
# buffer_width_um = 80 -> 10 bins (80 / 8)
# min_cluster_area_um2 = 1280 -> 20 bins
spe <- detectEdgeArtifacts_VisiumHD(spe, resolution = "8um")

# 16um data with defaults  
# buffer_width_um = 80 -> 5 bins (80 / 16)
# min_cluster_area_um2 = 1280 -> 5 bins
spe <- detectEdgeArtifacts_VisiumHD(spe, resolution = "16um")

# Custom parameters (physical units)
spe <- detectEdgeArtifacts_VisiumHD(
  spe, 
  resolution = "16um",
  buffer_width_um = 100,      # 100 $u$m buffer
  min_cluster_area_um2 = 2000 # 2000 $u$m^2 minimum 
)
} # }
```
