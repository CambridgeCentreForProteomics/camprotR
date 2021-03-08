context('combine_runs')

#### Setup ---------------------------------------------------------------------


#### Tests ---------------------------------------------------------------------
test_that("get_parsimony_pep2prot works on some SILAC peptide data", {
  expect_equal_to_reference(
    get_parsimony_pep2prot(list(
      system.file("testdata", "small_pep_silac.txt", package = "camprotR"),
      system.file("testdata", "small_pep_silac2.txt", package = "camprotR"))),
    "reference/get_parsimony_pep2prot.rds"
  )
})

test_that("single protein, matched", {
  expect_true(camprotR:::match_proteins('prot1', 'prot1'))
})

test_that("single protein, unmatched", {
  expect_false(camprotR:::match_proteins('prot1', 'prot2'))
})

test_that("multiple proteins, matched, same order", {
  expect_true(camprotR:::match_proteins('prot1; prot2', 'prot1; prot2'))
})

test_that("multiple proteins, matched, different order", {
  expect_true(camprotR:::match_proteins('prot1; prot2', 'prot2; prot1'))
})

test_that("multiple proteins, matched, different sep ", {
  expect_true(camprotR:::match_proteins('prot1:prot2', 'prot2; prot1',
                                        sep1=':', sep2='; '))
})

test_that("multiple proteins, unmatched", {
  expect_false(camprotR:::match_proteins('prot1; prot2', 'prot2; prot3'))
})

#### Sanity checks -------------------------------------------------------------
