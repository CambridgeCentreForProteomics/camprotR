# Load libraries
library(dplyr)

# All files from https://github.com/CambridgeCentreForProteomics/RBPMap_RNase/tree/master/pilot_3/raw
# Read in PSM data
'./data-raw/RBPMap_RNase_pilot_3/LFQ_OOPS_RNase_MetChloroform-(1)_PeptideGroups.txt'
lfq_oops_rnase <- read.delim(
  here::here("data-raw/RBPMap_RNase_pilot_3/LFQ_OOPS_RNase_MetChloroform-(1)_PeptideGroups.txt")
)

# Adding an artifical protein accessions column as this is needed for camprotR::parse_features
lfq_oops_rnase$Protein.Accessions <- lfq_oops_rnase$Master.Protein.Accessions

# Output .rda file
usethis::use_data(lfq_oops_rnase)
