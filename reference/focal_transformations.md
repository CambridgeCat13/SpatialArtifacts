# Apply morphological transformations to identify connected outlier regions

Performs a series of focal operations to clean and connect outlier
regions in spatial transcriptomics data through morphological
operations.

## Usage

``` r
focal_transformations(raster_object, min_cluster_size = 40)
```

## Arguments

- raster_object:

  A terra SpatRaster object with binary values where 1 indicates outlier
  spots and 0/NA indicates normal spots

- min_cluster_size:

  Numeric value specifying the minimum size for isolated clusters.
  Clusters smaller than this threshold will be filled (default: 40)

## Value

A processed SpatRaster object with cleaned and connected outlier regions

## Details

The function applies four sequential morphological operations:

1.  `3x3` fill: Fills spots completely surrounded by outliers

2.  `5x5` outline: Fills spots outlined by outliers in larger window

3.  Star pattern: Fills spots with outliers in cardinal directions

4.  Small cluster removal: Removes isolated normal regions below
    threshold

## Examples

``` r
library(terra)
#> terra 1.8.93
#> 
#> Attaching package: ‘terra’
#> The following objects are masked from ‘package:SummarizedExperiment’:
#> 
#>     distance, nearest, shift, trim, values, values<-, width
#> The following objects are masked from ‘package:GenomicRanges’:
#> 
#>     distance, gaps, nearest, shift, trim, values, values<-, width
#> The following objects are masked from ‘package:IRanges’:
#> 
#>     distance, gaps, nearest, shift, trim, width
#> The following objects are masked from ‘package:S4Vectors’:
#> 
#>     values, values<-, width
#> The following object is masked from ‘package:BiocGenerics’:
#> 
#>     width
#> The following object is masked from ‘package:generics’:
#> 
#>     interpolate
# Create a 5x5 mock raster object with an outlier (1)
m <- matrix(0, nrow = 5, ncol = 5)
m[2, 2] <- 1 # A single outlier
r <- terra::rast(m) # Use terra instead of raster

# Apply morphological cleaning
r_cleaned <- focal_transformations(r, min_cluster_size = 3)

# See the original and cleaned values
print(terra::values(r))
#>       lyr.1
#>  [1,]     0
#>  [2,]     0
#>  [3,]     0
#>  [4,]     0
#>  [5,]     0
#>  [6,]     0
#>  [7,]     1
#>  [8,]     0
#>  [9,]     0
#> [10,]     0
#> [11,]     0
#> [12,]     0
#> [13,]     0
#> [14,]     0
#> [15,]     0
#> [16,]     0
#> [17,]     0
#> [18,]     0
#> [19,]     0
#> [20,]     0
#> [21,]     0
#> [22,]     0
#> [23,]     0
#> [24,]     0
#> [25,]     0
print(terra::values(r_cleaned))
#>       lyr.1
#>  [1,]     0
#>  [2,]     0
#>  [3,]     0
#>  [4,]     0
#>  [5,]     0
#>  [6,]     0
#>  [7,]     1
#>  [8,]     0
#>  [9,]     0
#> [10,]     0
#> [11,]     0
#> [12,]     0
#> [13,]     0
#> [14,]     0
#> [15,]     0
#> [16,]     0
#> [17,]     0
#> [18,]     0
#> [19,]     0
#> [20,]     0
#> [21,]     0
#> [22,]     0
#> [23,]     0
#> [24,]     0
#> [25,]     0
```
