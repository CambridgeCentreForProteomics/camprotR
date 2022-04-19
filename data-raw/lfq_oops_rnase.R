# Load libraries
library(dplyr)

# All files from https://github.com/CambridgeCentreForProteomics/RBPMap_RNase/tree/master/pilot_3/raw/reexport
# Read in PSM data
lfq_oops_rnase_pep <- read.delim(
  here::here("data-raw/RBPMap_RNase_pilot_3/LFQ_OOPS_RNase_MetChloroform-(1)_PeptideGroups.txt")
)

random_proteins <- sample(unique(lfq_oops_rnase_pep$Master.Protein.Accessions), 100)
lfq_oops_rnase_pep <- lfq_oops_rnase_pep %>% filter(Master.Protein.Accessions %in% random_proteins)


# Output .rda file
usethis::use_data(lfq_oops_rnase_pep, overwrite=TRUE)

lfq_oops_rnase_psm <- read.delim(
  here::here("data-raw/RBPMap_RNase_pilot_3/LFQ_OOPS_RNase_MetChloroform-(1)_PSMs.txt")
)


lfq_oops_rnase_psm <- lfq_oops_rnase_psm %>% filter(Master.Protein.Accessions %in% random_proteins)

# Output .rda file
usethis::use_data(lfq_oops_rnase_psm, overwrite=TRUE)
