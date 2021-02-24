# Read in SILAC PSM level data
silac <- read.delim(
  here::here("data-raw/H1_PSMs.txt")
)

# Set seed for reproducible sampling
set.seed(2021)

# Randomly sample 1000 rows
silac_filt <- silac[sample(nrow(silac), 1000), ]

# Select the columns we need
small_psm_silac <- silac_filt %>%
  select(
    Annotated.Sequence,
    Modifications,
    Isolation.Interference.in.Percent,
    Quan.Channel,
    Precursor.Abundance
  )

# Output .rda file
save(small_psm_silac, compress = "bzip2", file = "inst/testdata/small_psm_silac.rda")
