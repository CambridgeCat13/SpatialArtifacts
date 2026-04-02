# inst/scripts/make-data.R
#
# Script to reproduce the spe_vignette dataset included in SpatialArtifacts.
#
# Data source:
#   Human hippocampus Visium data (sample V11L05-335_C1)
#   from Thompson et al. (2025), Nature Neuroscience.
#   https://doi.org/10.1038/s41593-025-02022-0
#
#   Original Space Ranger output located at:
#   /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/
#     processed-data/01_spaceranger/spaceranger-all/V11L05-335_C1/outs/
#
# Note: Raw data is not publicly available via this internal path.
#   Users may access the processed data via the humanHippocampus2024
#   Bioconductor ExperimentHub package:
#   https://bioconductor.org/packages/humanHippocampus2024
#   Or contact the Hicks Lab at Johns Hopkins University.

library(SpatialExperiment)
library(scuttle)
library(scran)
library(Matrix)
library(usethis)

# Load raw Space Ranger output
sample_path <- "/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/01_spaceranger/spaceranger-all/V11L05-335_C1/outs/"

spe <- read10xVisium(
  samples = sample_path,
  sample_id = "V11L05-335_C1",
  type = "HDF5",
  data = "filtered"
)

# 2. Compute QC metrics (used by detectEdgeArtifacts via 'sum' column)
#    This must be done BEFORE gene subsetting so metrics reflect all genes
spe <- addPerCellQC(spe)

# 3. Subset to informative genes to reduce package size
#    Filter zero-variance genes, retaining ~12,971 informative genes
spe_norm <- logNormCounts(spe)
gene_var <- modelGeneVar(spe_norm)
top_genes <- getTopHVGs(gene_var, var.threshold = 0)
spe <- spe[top_genes, ]

# Remove image data to meet Bioconductor 5MB size limit
imgData(spe) <- NULL

# Sparsify counts matrix
assay(spe, "counts") <- Matrix::Matrix(
  as.matrix(assay(spe, "counts")), sparse = TRUE
)

spe_vignette <- spe
usethis::use_data(spe_vignette, compress = "xz", overwrite = TRUE)