context('PSM to peptide style conversion')

test_that("N-Term, multiople Carbamidomethyl", {
  expect_equal(psm_to_peptide_style_modifications(
    "N-Term(Prot)(Acetyl); C2(Carbamidomethyl); C16(Carbamidomethyl)"),
    "1xAcetyl [N-Term]; 2xCarbamidomethyl [C2; C16]")
})

