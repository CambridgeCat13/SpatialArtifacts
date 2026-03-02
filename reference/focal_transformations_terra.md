# Morphological transformations using terra (optimized for VisiumHD)

Terra-only implementation of focal operations for connecting outlier
regions.

## Usage

``` r
focal_transformations_terra(r, min_cluster_size = 5)
```

## Arguments

- r:

  A
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object with binary values (1=outlier, 0=normal).

- min_cluster_size:

  Minimum size (in bins) for small hole removal.

## Value

A
[`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
object after applying focal transformations (fill, outline, star) and
small hole removal.

## Examples

``` r
# Create a dummy binary SpatRaster for demonstration
if (requireNamespace("terra", quietly = TRUE)) {
  r <- terra::rast(matrix(c(0, 0, 0, 0, 1, 0, 0, 0, 0), nrow = 3),
    extent = terra::ext(0, 3, 0, 3)
  )
  # Run the focal transformations
  # result <- focal_transformations_terra(r, min_cluster_size = 5)
}
```
