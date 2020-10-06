context("Parsing PSM level PD output for TMT data")

#### Setup ---------------------------------------------------------------------
psm_data <- read.delim(
  system.file("extdata", "TMT_PSMs.txt.gz", package = "camprotR"),
  stringsAsFactors=FALSE)

crap_fasta_inf <- system.file("extdata", "cRAP_FullIdentifiers.fasta.gz",
                              package = "camprotR")

# Load the cRAP FASTA used for the PD search
crap.fasta <- Biostrings::fasta.index(crap_fasta_inf, seqtype = "AA")

# Extract the non cRAP UniProt accessions associated with each cRAP protein
crap.accessions <- unlist(
  stringr::str_extract_all(crap.fasta$desc, "(?<=\\|).*?(?=\\|)"))

#### Tests ---------------------------------------------------------------------
out <- camprotR::parse_features(psm_data, TMT=TRUE, level='PSM',
                                crap_proteins=crap.accessions)
expect_equal_to_reference(out, file='parsed_tmt_psm.rds')
