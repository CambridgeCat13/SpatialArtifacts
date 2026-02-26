# Identify all problem areas in spatial transcriptomics data

Detects and characterizes all clusters of low-quality spots (problem
areas) in tissue sections, including both edge and interior artifacts.

## Usage

``` r
problemAreas(
  .xyz,
  offTissue,
  uniqueIdentifier = NA,
  shifted = FALSE,
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

- uniqueIdentifier:

  Character string used as prefix for cluster IDs. If NA, "X" will be
  used (default: NA)

- shifted:

  Logical indicating whether to apply coordinate adjustment for
  hexagonal arrays (default: FALSE)

- min_cluster_size:

  Numeric value for minimum cluster size in morphological operations
  (default: 40)

## Value

Data frame with the following columns:

- spotcode: Spot identifier

- clumpID: Unique cluster identifier (format: "prefix_number")

- clumpSize: Number of spots in the cluster

Returns empty data frame if no problem areas are found.

## Details

This function identifies ALL connected components of outlier spots, not
just those touching edges. Each cluster is assigned a unique ID and its
size is calculated. This enables downstream filtering based on cluster
characteristics.

## Examples

``` r
# 1. Create a 5x5 grid of mock spot data
spot_data <- data.frame(
  array_row = rep(1:5, each = 5),
  array_col = rep(1:5, times = 5),
  outlier = FALSE
)
rownames(spot_data) <- paste0("spot", 1:25)

# 2. Create an "artifact" by flagging a 2x2 area as outliers
spot_data$outlier[spot_data$array_row %in% 2:3 & spot_data$array_col %in% 2:3] <- TRUE

# 3. Define off-tissue spots
offTissue <- character(0)

# 4. Identify all problem areas
problem_df <- problemAreas(
  spot_data,
  offTissue = offTissue,
  uniqueIdentifier = "Sample1",
  min_cluster_size = 3
)
print(problem_df)
#>   spotcode   clumpID clumpSize
#> 1    spot8 Sample1_1         4
#> 2    spot7 Sample1_1         4
#> 3   spot13 Sample1_1         4
#> 4   spot12 Sample1_1         4
```
