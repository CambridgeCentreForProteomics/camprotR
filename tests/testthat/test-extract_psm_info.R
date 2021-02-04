context('extract_psm_info')

#### Setup ---------------------------------------------------------------------

#### Tests ---------------------------------------------------------------------
test_that("Extracting sequenced and interference data from PSM", {
  data('silac_psm_data', package='camprotR')
  out <-  silac_psm_seq_int(silac_psm_data, sequence_col='Annotated.Sequence') %>%
    dplyr::arrange(Annotated.Sequence, Modifications)

  #expect_equal_to_reference(
  #  dplyr::arrange(output, Annotated.Sequence, Modifications), # sort order appears to be env-specific
  #  file='psm_sequenced.rds')
  ref <- readRDS('psm_sequenced.rds') %>%
    dplyr::arrange(Annotated.Sequence, Modifications)
  expect_equal(out, ref)
})

#### Sanity checks -------------------------------------------------------------
