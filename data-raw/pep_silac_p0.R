library(dplyr)

# Read in SILAC passage 0 peptide data
# (i.e. cells just in R0K0 media)
pepg <- read.delim(
  here::here("data-raw/Molm_13_P0_PeptideGroups.txt")
)

# Set seed for reproducible sampling
set.seed(2021)

# Get all peptides for a random sample of 100 Master proteins
pep_silac_p0 <- pepg %>%
  filter(Master.Protein.Accessions %in% sample(unique(Master.Protein.Accessions), 100))

# Output .rda file
usethis::use_data(pep_silac_p0)
