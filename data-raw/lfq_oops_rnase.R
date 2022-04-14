# Load libraries
library(dplyr)

# All files from https://github.com/CambridgeCentreForProteomics/RBPMap_RNase/tree/master/pilot_3/raw/reexport
# Read in PSM data
lfq_oops_rnase_pep <- read.delim(
  here::here("data-raw/RBPMap_RNase_pilot_3/LFQ_OOPS_RNase_MetChloroform-(1)_PeptideGroups.txt")
)

# Output .rda file
usethis::use_data(lfq_oops_rnase_pep)

lfq_oops_rnase_psm <- read.delim(
  here::here("data-raw/RBPMap_RNase_pilot_3/LFQ_OOPS_RNase_MetChloroform-(1)_PSMs.txt")
)

# Output .rda file
usethis::use_data(lfq_oops_rnase_psm)
