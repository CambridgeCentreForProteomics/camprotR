# Load libraries
library(dplyr) # I used version 1.0.4
library(MSnbase) # I used version 2.16.1

# Read in PSM data
tmt_total <- read.delim(
  here::here("data-raw/anja_lopit_total_proteome_rep1_PSMs.txt")
)

# Filter for first 100 PSMs
tmt_total_filt <- tmt_total[1:100, ]

# Extract expression data (i.e. PSM abundances in each TMT channel)
abundance_cols <- colnames(tmt_total_filt)[grepl('Abundance.', colnames(tmt_total_filt))]
e_data <- as.matrix(tmt_total_filt[, abundance_cols])

# Extract feature data that we want (i.e. metadata about each PSM)
f_data <- tmt_total_filt[, c("Master.Protein.Accessions",
                               "Isolation.Interference.in.Percent",
                               "Average.Reporter.SN")]

# Construct MSnSet
small_psm_tmt_total <- MSnbase::MSnSet(exprs = output_e, fData = output_f)

# Output .rda file
save(small_psm_tmt_total, compress = "bzip2", file = here::here("tests/testthat/data-test/small_psm_tmt_total.rda"))
