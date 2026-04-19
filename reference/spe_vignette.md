# Example SpatialExperiment for Vignettes

A lightweight
[`SpatialExperiment`](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment.html)
object derived from a human hippocampus Visium sample, used for
demonstrating the `SpatialArtifacts` artifact detection workflow.

## Format

A
[`SpatialExperiment`](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment.html)
object with 12,971 genes and 4,965 spots, containing one sample
(V11L05-335_C1). The `colData` includes precomputed QC metrics (e.g.,
`sum`, `detected`) from
[`addPerCellQC`](https://rdrr.io/pkg/scuttle/man/addPerCellQCMetrics.html).

## Source

Derived from human hippocampus Visium data (sample V11L05-335_C1) from
Thompson et al. (2025). The raw Space Ranger output was accessed
internally via the spatialHPC project (LIBD4035) on the JHPCE cluster.
The same dataset is publicly available via the `humanHippocampus2024`
Bioconductor ExperimentHub package
(<https://bioconductor.org/packages/humanHippocampus2024>). A script to
reproduce this object is provided in `inst/scripts/make-data.R`.

## Value

A
[`SpatialExperiment`](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment.html)
object.

## Details

The object was derived from human hippocampus Visium data (sample
V11L05-335_C1) from the spatialHPC project (LIBD4035). To meet
Bioconductor package size requirements (\<5 MB), the object was subset
to highly variable genes using `scran::getTopHVGs()`, image data was
removed, and the counts matrix was stored as a sparse matrix with XZ
compression. QC metrics in `colData` were computed prior to gene
subsetting and remain accurate for the full spot set.

## References

Thompson, J.R., Nelson, E.D., Tippani, M. et al. (2025). An integrated
single-nucleus and spatial transcriptomics atlas reveals the molecular
landscape of the human hippocampus. *Nature Neuroscience* 28, 1990–2004.
[doi:10.1038/s41593-025-02022-0](https://doi.org/10.1038/s41593-025-02022-0)

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
