test_that("classifyEdgeArtifacts returns expected columns", {
  data('spe_vignette', package = 'SpatialArtifacts')
  SummarizedExperiment::assay(spe_vignette, 'counts') <- as.matrix(SummarizedExperiment::assay(spe_vignette, 'counts'))
  SummarizedExperiment::colData(spe_vignette)$sum_umi <- SummarizedExperiment::colData(spe_vignette)$sum

  # Run detect first
  spe_detected <- detectEdgeArtifacts(
    spe_vignette,
    platform = 'visium',
    qc_metric = 'sum_umi',
    batch_var = 'sample_id',
    shifted = FALSE
  )

  #  Run classify
  spe_classified <- classifyEdgeArtifacts(
    spe_detected,
    qc_metric = 'sum_umi',
    min_spots = 20,
    name = 'edge_artifact'
  )

  # 5. Validate
  expect_s4_class(spe_classified, 'SpatialExperiment')
  expect_true('edge_artifact_classification' %in% 
                colnames(SummarizedExperiment::colData(spe_classified)))
  
  valid_classes <- c('not_artifact', 'large_edge_artifact', 'small_edge_artifact',
                     'large_interior_artifact', 'small_interior_artifact')
  expect_true(all(spe_classified$edge_artifact_classification %in% valid_classes))
  expect_true('not_artifact' %in% spe_classified$edge_artifact_classification)
})

test_that("classifyEdgeArtifacts errors without required columns", {
  data('spe_vignette', package = 'SpatialArtifacts')
  
  expect_error(
    classifyEdgeArtifacts(spe_vignette, name = 'edge_artifact'),
    "Missing required columns"
  )
})