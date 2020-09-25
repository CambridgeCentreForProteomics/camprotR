context('Extracting sequenced and interference data from PSM')

test_that("Extract info", {
  data('silac_psm_data', package='camprotR')
  output <-  silac_psm_seq_int(silac_psm_data, sequence_col='Annotated.Sequence')

  expect_equal_to_reference(
    dplyr::arrange(output, Annotated.Sequence, Modifications), # sort order appears to be env-specific
    file='psm_sequenced.rds')
})
