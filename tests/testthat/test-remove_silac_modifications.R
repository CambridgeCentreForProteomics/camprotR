context('Remove SILAC modifications')

test_that("PSM", {
  expect_equal(remove_silac_modifications(
    "N-Term(Prot)(Acetyl); R15(Label:13C(6)15N(4))", level='psm'),
    "N-Term(Prot)(Acetyl)")
})


test_that("Peptide", {
  expect_equal(remove_silac_modifications(
    "1xCarbamidomethyl [C9]; 1xLabel:13C(6)15N(2) [K12]", level='peptide'),
    "1xCarbamidomethyl [C9]")
})
