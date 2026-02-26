# Detect edge dryspots using spatial clustering

Identifies clusters of low-quality spots that touch tissue boundaries,
indicating potential edge dryspot artifacts from incomplete reagent
coverage.

## Usage

``` r
clumpEdges(
  .xyz,
  offTissue,
  shifted = FALSE,
  edge_threshold = 0.75,
  min_cluster_size = 40
)
```

## Arguments

- .xyz:

  Data frame with spot coordinates and QC metrics. Must contain:

  - Column 1: array_row coordinates

  - Column 2: array_col coordinates

  - Column 3: Binary outlier indicator (TRUE/FALSE or 1/0)

  - Row names: Spot identifiers

- offTissue:

  Character vector of spot identifiers that are off-tissue

- shifted:

  Logical indicating whether to apply coordinate adjustment for
  hexagonal arrays (default: FALSE)

- edge_threshold:

  Numeric value between 0 and 1 specifying the minimum proportion of
  border coverage required for edge detection (default: 0.75)

- min_cluster_size:

  Numeric value for minimum cluster size in morphological operations
  (default: 40)

## Value

Character vector of spot identifiers classified as edge dryspots

## Details

The function performs the following steps:

1.  Converts spot coordinates to raster format

2.  Applies morphological transformations to connect outlier regions

3.  Identifies connected components (clusters)

4.  Checks each tissue border (N, S, E, W) for cluster coverage

5.  Returns spots from clusters that exceed edge_threshold on any border

## Examples

``` r
# 1. Create a 5x5 grid of mock spot data
spot_data <- data.frame(
  array_row = rep(1:5, each = 5),
  array_col = rep(1:5, times = 5),
  outlier = FALSE
)
rownames(spot_data) <- paste0("spot", 1:25)

# 2. Create an "edge artifact" by flagging the first row as outliers
spot_data$outlier[spot_data$array_row == 1] <- TRUE

# 3. Define off-tissue spots (none in this simple case)
offTissue <- character(0)

# 4. Detect edge dryspots
edge_spots <- clumpEdges(
  spot_data,
  offTissue = offTissue,
  edge_threshold = 0.5,
  min_cluster_size = 3
)
print(edge_spots)
#> [1] "spot1" "spot2" "spot3" "spot4" "spot5"
```
