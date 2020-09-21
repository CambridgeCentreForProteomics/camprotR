context('Extracting sequenced and interference data from PSM')

test_that("Extract info", {
  data('silac_psm_data', package='camprotR')
  expect_equal_to_reference(
    arrange(silac_psm_seq_int(silac_psm_data, sequence_col='Annotated.Sequence'), Annotated.Sequence),
    file='psm_sequenced.rds')
})
