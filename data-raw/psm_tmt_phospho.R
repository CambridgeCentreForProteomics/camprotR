# Read in PSM data
tmt_phospho <- read.delim(
  here::here("data-raw/Anja_phospho_lopit_phosphoproteome_rep2_PSMs.txt")
)

# Set seed for reproducible sampling
set.seed(2021)

# Randomly sample 5000 rows
psm_tmt_phospho <- tmt_phospho[sample(nrow(tmt_phospho), 5000), ]

# Output .rda file
usethis::use_data(psm_tmt_phospho)
