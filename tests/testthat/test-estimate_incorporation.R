test_that("estimate_incorporation() works", {
  # estimate incorporation and output to incorporation folder
  out <- estimate_incorporation(
    psm_input = psm_silac_p4,
    peptide_input = pep_silac_p4,
    crap_fasta = system.file("extdata/cRAP_20190401.fasta.gz", package = "camprotR"),
    mix = 1
  )

  # check plot outputs
  vdiffr::expect_doppelganger(
    "estimate-incorporation-HL-correlation",
    out$HL_correlation
  )

  vdiffr::expect_doppelganger(
    "estimate-incorporation-peptide",
    out$peptide_incorporation
  )

  vdiffr::expect_doppelganger(
    "estimate-incorporation-protein",
    out$protein_incorporation
  )

  # check table output
  withr::with_tempfile("tf", {
    # save output to tempfile
    write.table(out$incorporation_table, file = tf,
                sep = "\t", row.names = FALSE, col.names = TRUE)

    expect_snapshot_file(tf, "estimate-incorporation-table")
  })
})
