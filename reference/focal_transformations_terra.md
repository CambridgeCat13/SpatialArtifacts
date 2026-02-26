# Morphological transformations using terra (optimized for VisiumHD)

Terra-only implementation of focal operations for connecting outlier
regions.

## Usage

``` r
focal_transformations_terra(r, min_cluster_size = 5)
```

## Arguments

- r:

  Terra SpatRaster with binary values (1=outlier, 0=normal)

- min_cluster_size:

  Minimum size for small hole removal

## Value

Processed SpatRaster
