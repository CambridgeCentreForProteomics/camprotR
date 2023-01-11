#### Tests ---------------------------------------------------------------------
test_that("Modification style conversion works, N-Term, multiple Carbamidomethyl", {
  expect_equal(psm_to_peptide_style_modifications(
    "N-Term(Prot)(Acetyl); C2(Carbamidomethyl); C16(Carbamidomethyl)"),
    "1xAcetyl [N-Term]; 2xCarbamidomethyl [C2; C16]")
})

test_that("PSM style SILAC modifications can be removed", {
  expect_equal(remove_silac_modifications(
    "N-Term(Prot)(Acetyl); R15(Label:13C(6)15N(4))", level='psm'),
    "N-Term(Prot)(Acetyl)")
})


test_that("Peptide style SILAC modifications can be removed", {
  expect_equal(remove_silac_modifications(
    "1xCarbamidomethyl [C9]; 1xLabel:13C(6)15N(2) [K12]", level='peptide'),
    "1xCarbamidomethyl [C9]")
})

test_that("get_psm_silac_mod_regex works with no args", {
  df <- get_psm_silac_mod_regex()

  expect_true(inherits(df, "data.frame"))
  expect_equal(colnames(df), c("name", "desc", "regex"))
})

test_that("get_psm_silac_mod_regex works with args", {
  regex <- get_psm_silac_mod_regex("K_2H4")

  expect_equal(regex, 'K\\d{1,2}\\(Label:2H\\(4\\)\\)')
})

#### Sanity checks -------------------------------------------------------------
test_that("remove_silac_modifications() errors if level is nonsense", {
  expect_error(remove_silac_modifications(
    "1xCarbamidomethyl [C9]; 1xLabel:13C(6)15N(2) [K12]", level='banana'),
    "level must be psm or peptide")
})
