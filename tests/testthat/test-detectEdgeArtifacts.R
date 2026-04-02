test_that('detectEdgeArtifacts returns expected columns', {
  # 1. Load test data
  data('spe_vignette', package = 'SpatialArtifacts')
  
  # 2. Prepare data
  SummarizedExperiment::assay(spe_vignette, 'counts') <- as.matrix(SummarizedExperiment::assay(spe_vignette, 'counts'))
  SummarizedExperiment::colData(spe_vignette)$sum_umi <- SummarizedExperiment::colData(spe_vignette)$sum
  
  # 3. Run function
  spe_res <- detectEdgeArtifacts(
    spe_vignette, 
    platform = 'visium', 
    qc_metric = 'sum_umi', 
    batch_var = 'sample_id', 
    shifted = FALSE
  )
  
  # 4. Validate
  expect_s4_class(spe_res, 'SpatialExperiment')
  expect_true('edge_artifact_edge' %in% colnames(SummarizedExperiment::colData(spe_res)))
  expect_true('edge_artifact_problem_id' %in% colnames(SummarizedExperiment::colData(spe_res)))
  expect_true('edge_artifact_problem_size' %in% colnames(SummarizedExperiment::colData(spe_res)))
})
