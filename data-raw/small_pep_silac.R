# Read in SILAC peptide level data
silac <- read.delim(
  here::here("data-raw/H1_PeptideGroups.txt")
)

# Set seed for reproducible sampling
set.seed(2021)

# Randomly sample 1000 rows
silac_filt <- silac[sample(nrow(silac), 1000), ]

# Select the columns we need
output <- silac_filt %>%
  select(
    Sequence,
    Modifications,
    Master.Protein.Accessions,
    Protein.Accessions,
    Number.of.Protein.Groups,
    Quan.Info,
    Abundances.Grouped.Light,
    Abundances.Grouped.Heavy
  )

# Output text file
readr::write_delim(output, "inst/testdata/small_pep_silac.txt", delim = "\t")
