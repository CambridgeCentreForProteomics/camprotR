context('Extracting sequenced and interference data from PSM')

test_that("Extract info", {
  data('silac_psm_data', package='camprotR')
  output <-  silac_psm_seq_int(silac_psm_data, sequence_col='Annotated.Sequence')

  ref_file <- 'psm_sequenced.rds'
  print(file.exists(ref_file))
  ref <- readRDS(ref_file)
  print(nrow(ref))
  print(nrow(dplyr::arrange(output, Annotated.Sequence, Modifications)))
  print(ref[1,])
  print(dplyr::arrange(output, Annotated.Sequence, Modifications)[1,])
  stop()

  expect_equal_to_reference(
    dplyr::arrange(output, Annotated.Sequence, Modifications), # sort order appears to be env-specific
    file='psm_sequenced.rds')
})
