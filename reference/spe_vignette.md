# Mock SpatialExperiment for Vignettes

A minimal
[`SpatialExperiment`](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment.html)
object used for demonstrating the spatial artifact detection workflow.

## Format

A
[`SpatialExperiment`](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment.html)
object.

## Source

Internal simulation

## Value

A
[`SpatialExperiment`](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment.html)
object.

## Examples

``` r
data(spe_vignette)
spe_vignette
#> class: SpatialExperiment 
#> dim: 12971 4965 
#> metadata(0):
#> assays(1): counts
#> rownames(12971): ENSG00000244734 ENSG00000123560 ... ENSG00000254089
#>   ENSG00000125900
#> rowData names(0):
#> colnames(4965): AAACAACGAATAGTTC-1 AAACAAGTATCTCCCA-1 ...
#>   TTGTTTGTATTACACG-1 TTGTTTGTGTAAATTC-1
#> colData names(6): sample_id in_tissue ... sum detected
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
#> imgData names(0):
```
