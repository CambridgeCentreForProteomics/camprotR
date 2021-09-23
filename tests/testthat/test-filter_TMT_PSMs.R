save_update_average_sn <- function(msnset) {
  path <- tempfile(tmpdir = test_path("testdata"), fileext = ".txt")
  out <- update_average_sn(
    obj = msnset,
    sn_col = "Average.Reporter.SN"
  )
  write.table(out, file = path,
              sep = "\t", row.names = FALSE, col.names = TRUE)

  return(path)
}

test_that("update_average_sn() works", {
  load(test_path("testdata/small_psm_tmt_total.rda"))

  expect_snapshot_file(
    save_update_average_sn(
      msnset = small_psm_tmt_total
    ),
    "update_average_sn.txt"
  )
})

save_filter_TMT_PSMs <- function(msnset, inter_thresh, sn_thresh, output = c("exprs", "fData")) {
  path <- tempfile(tmpdir = test_path("testdata"), fileext = ".txt")
  out <- filter_TMT_PSMs(
    obj = msnset,
    inter_thresh = inter_thresh,
    sn_thresh = sn_thresh,
    master_protein_col = "Master.Protein.Accessions",
    inter_col = "Isolation.Interference.in.Percent",
    sn_col = "Average.Reporter.SN",
    verbose = TRUE
  )

  if (output == "exprs") {
    write.table(MSnbase::exprs(out), file = path, sep = "\t",
                row.names = FALSE, col.names = TRUE)
  } else if (output == "fData") {
    write.table(MSnbase::fData(out), file = path, sep = "\t",
                row.names = FALSE, col.names = TRUE)
  }

  return(path)
}

test_that("filter_TMT_PSMs() works", {
  load(test_path("testdata/small_psm_tmt_total.rda"))

  expect_snapshot_file(
    save_filter_TMT_PSMs(
      msnset = small_psm_tmt_total,
      inter_thresh = 0,
      sn_thresh = 10,
      output = "exprs"
    ),
    "filter_TMT_PSMs_exprs.txt"
  )

  expect_snapshot_file(
    save_filter_TMT_PSMs(
      msnset = small_psm_tmt_total,
      inter_thresh = 0,
      sn_thresh = 10,
      output = "fData"
    ),
    "filter_TMT_PSMs_fData.txt"
  )
})
