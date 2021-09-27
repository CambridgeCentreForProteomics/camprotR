snap_update_average_sn <- function(msnset) {
  # make a tempfile to compare against snapshot file
  path <- tempfile(tmpdir = test_path("testdata"), fileext = ".txt")

  # run function to test and capture output
  out <- update_average_sn(
    obj = msnset,
    sn_col = "Average.Reporter.SN"
  )

  # save output to tempfile
  write.table(out, file = path,
              sep = "\t", row.names = FALSE, col.names = TRUE)

  # return path of tempfile
  return(path)
}

test_that("update_average_sn() works", {
  # load in test MSnSet
  load(test_path("testdata/small_psm_tmt_total.rda"))

  # compare tempfile to snapshot file
  expect_snapshot_file(
    snap_update_average_sn(
      msnset = small_psm_tmt_total
    ),
    "update_average_sn.txt"
  )
})

snap_filter_TMT_PSMs <- function(msnset, inter_thresh, sn_thresh, output = c("exprs", "fData")) {
  # make a tempfile to compare against snapshot file
  path <- tempfile(tmpdir = test_path("testdata"), fileext = ".txt")

  # run function to test and capture output
  out <- filter_TMT_PSMs(
    obj = msnset,
    inter_thresh = inter_thresh,
    sn_thresh = sn_thresh,
    master_protein_col = "Master.Protein.Accessions",
    inter_col = "Isolation.Interference.in.Percent",
    sn_col = "Average.Reporter.SN",
    verbose = TRUE
  )

  # extract data from output MSnSet and save to tempfile
  if (output == "exprs") {
    write.table(MSnbase::exprs(out), file = path, sep = "\t",
                row.names = FALSE, col.names = TRUE)
  } else if (output == "fData") {
    write.table(MSnbase::fData(out), file = path, sep = "\t",
                row.names = FALSE, col.names = TRUE)
  }

  # return path of tempfile
  return(path)
}

test_that("filter_TMT_PSMs() works", {
  # load in test MSnSet
  load(test_path("testdata/small_psm_tmt_total.rda"))

  # compare tempfile (expression data) to snapshot file
  expect_snapshot_file(
    snap_filter_TMT_PSMs(
      msnset = small_psm_tmt_total,
      inter_thresh = 0,
      sn_thresh = 10,
      output = "exprs"
    ),
    "filter_TMT_PSMs_exprs.txt"
  )

  # compare tempfile (feature data) to snapshot file
  expect_snapshot_file(
    snap_filter_TMT_PSMs(
      msnset = small_psm_tmt_total,
      inter_thresh = 0,
      sn_thresh = 10,
      output = "fData"
    ),
    "filter_TMT_PSMs_fData.txt"
  )
})
