context('Adding peptide positions')

# set up #
input <- data.frame('protein'=c('L0R819', 'P98196', 'P14618'),
                    'Sequence'=c('GTRPEDSSVLIPTDNSTPHK',
                                 'QLWPTATER',
                                 'LNFSHGTHEYHAETIK'))


expected_output <- input
expected_output$peptide_start <- c(5, 1099, NA)
expected_output$peptide_end <- c(24, 1107, NA)

fasta <- system.file("extdata", "reference.fasta", package = "camprotR")
##########
# tests #
test_that("Simple input", {
  expect_equal(add_peptide_positions(input, fasta, master_protein_col = "protein"),
               expected_output)
})

