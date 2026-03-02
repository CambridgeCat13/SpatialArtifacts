# Detect problem areas in VisiumHD data using morphological processing

Terra-optimized implementation for VisiumHD. Uses morphological
operations to connect outlier regions and identify connected components.

## Usage

``` r
problemAreas_WithMorphology_terra(
  .xyz,
  uniqueIdentifier = NA,
  min_cluster_size = 5,
  resolution = "8um"
)
```

## Arguments

- .xyz:

  Data frame with columns: x, y, outlier (binary 0/1)

- uniqueIdentifier:

  Character string for cluster ID prefix

- min_cluster_size:

  Minimum cluster size in bins (default: 5)

- resolution:

  VisiumHD resolution: "8um" or "16um" (default: "8um")

## Value

Data frame with columns: spotcode, clumpID, clumpSize

## Examples

``` r
library(terra)

# 1. Create a mock VisiumHD-like coordinate dataframe
# 10x10 grid
coords <- expand.grid(x = 1:10, y = 1:10)
.xyz <- data.frame(
  x = coords$x,
  y = coords$y,
  outlier = 0
)
rownames(.xyz) <- paste0("bin_", 1:100)

# 2. Create a "problem area" (cluster of outliers)
# Make a 3x3 block of outliers in the center
center_mask <- .xyz$x >= 4 & .xyz$x <= 6 & .xyz$y >= 4 & .xyz$y <= 6
.xyz$outlier[center_mask] <- 1

# 3. Run detection
clusters <- problemAreas_WithMorphology_terra(
  .xyz,
  uniqueIdentifier = "TEST",
  min_cluster_size = 2,
  resolution = "16um"
)
#>     Clustering 9 outlier bins...
#>     Creating raster: 10 rows * 10 cols
#>     Applying morphological transformations...
#>     [16um Mode] Adjusted min_cluster_size: 2 -> 2
#>     Clustering outlier regions...
#>     Found 9 clustered cells
#>     Identified 1 unique clusters
#>     Kept 1 clusters (>=2 bins): 9 bins total

# Check results
print(clusters)
#>   spotcode clumpID clumpSize
#> 1   bin_34  TEST_1         9
#> 2   bin_44  TEST_1         9
#> 3   bin_54  TEST_1         9
#> 4   bin_35  TEST_1         9
#> 5   bin_45  TEST_1         9
#> 6   bin_55  TEST_1         9
#> 7   bin_36  TEST_1         9
#> 8   bin_46  TEST_1         9
#> 9   bin_56  TEST_1         9
```
