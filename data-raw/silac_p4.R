library(dplyr)

# Read in SILAC passage 4 peptide data
# 50:50 mix of light:heavy cells (R0K0:R10K8)
pepg <- read.delim(
  here::here("data-raw/Molm_13_P4_PeptideGroups.txt")
)

# Read in corresponding PSM data
silac_p4 <- read.delim(
  here::here("data-raw/Molm_13_P4_PSMs.txt")
)

# Remove extraneous periods from colnames
colnames(pepg) <- camprotR::remove_dots(colnames(pepg))
colnames(silac_p4) <- camprotR::remove_dots(colnames(silac_p4))

# Set seed for reproducible sampling
set.seed(2021)

# Get a random sample of 100 Master proteins
prot_sample <- sample(unique(pepg$Master.Protein.Accessions), 100)

# Get all peptides for this random sample
pep_silac_p4 <- filter(pepg, Master.Protein.Accessions %in% prot_sample)

# Get all PSMs for this random sample
psm_silac_p4 <- filter(silac_p4, Master.Protein.Accessions %in% prot_sample)

# Output peptide .rda file
usethis::use_data(pep_silac_p4, overwrite = TRUE)

# Output psm .rda file
usethis::use_data(psm_silac_p4, overwrite = TRUE)
