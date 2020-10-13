context('Adding site sequence')

# set up #
input <- data.frame('protein'=c('L0R819', 'P98196', 'P14618'),
                    'ptm_position'=c('20', '100', '40'))


expected_output <- input
expected_output$site_seq <- c('VLIPTDNsTPHKEDL', 'PVTSGLPlFFVITVT', NA)

fasta <- system.file("extdata", "reference.fasta", package = "camprotR")
##########
# tests #
test_that("Simple input", {
  expect_equal(add_site_sequence(input, fasta, master_protein_col='protein'),
               expected_output)
})
