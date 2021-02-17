context("parse_features")

#### Setup ---------------------------------------------------------------------
psm_data <- psm_tmt_total

crap_fasta_inf <- system.file("extdata", "cRAP_FullIdentifiers.fasta.gz",
                              package = "camprotR")

# Load the cRAP FASTA used for the PD search
crap.fasta <- Biostrings::fasta.index(crap_fasta_inf, seqtype = "AA")

# Extract the non cRAP UniProt accessions associated with each cRAP protein
crap.accessions <- regmatches(
  crap.fasta$desc,
  gregexpr("(?<=\\|).*?(?=\\|)", crap.fasta$desc, perl = TRUE)
) %>% unlist()

#### Tests ---------------------------------------------------------------------
test_that("parse_features works", {
  expect_equal_to_reference(
    parse_features(psm_data, TMT=TRUE, level='PSM',
                   crap_proteins=crap.accessions),
    file='reference/parsed_psm_tmt_total.rds'
  )
})

#### Sanity checks -------------------------------------------------------------
test_that("Error is produced if crap_proteins is not given", {
  expect_error(
    parse_features(psm_data, TMT = TRUE, level = "PSM"),
    "must supply the crap_proteins argument to filter cRAP proteins"
  )
})
