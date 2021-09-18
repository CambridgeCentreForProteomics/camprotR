library(dplyr)

# Read in SILAC passage 4 peptide data
# 50:50 mix of light:heavy cells (R0K0:R10K8)
pepg <- read.delim(
  here::here("data-raw/Molm_13_P4_PeptideGroups.txt")
)

# Set seed for reproducible sampling
set.seed(2021)

# Get all peptides for a random sample of 100 Master proteins
pep_silac_p4 <- pepg %>%
  filter(Master.Protein.Accessions %in% sample(unique(Master.Protein.Accessions), 100))

# Output .rda file
usethis::use_data(pep_silac_p4)
