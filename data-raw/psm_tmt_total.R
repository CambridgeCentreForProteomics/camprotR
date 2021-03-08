# Load libraries
library(dplyr)

# Read in PSM data
tmt_total <- read.delim(
  here::here("data-raw/anja_lopit_total_proteome_rep1_PSMs.txt")
)

# Filter for first 5000 PSMs
tmt_total_filt <- tmt_total[1:5000, ]

# Read in raw file metadata and extract needed info
tmt_total_meta <- read.delim(
  here::here("data-raw/anja_lopit_total_proteome_rep1_InputFiles.txt")
) %>%
  mutate(Spectrum.File = basename(gsub("\\\\", "/", File.Name))) %>%
  select(Study.File.ID, Spectrum.File)

# Process PSM table
psm_tmt_total <- tmt_total_filt %>%
  left_join(tmt_total_meta, by = c("File.ID" = "Study.File.ID")) %>%
  relocate(Spectrum.File, .before = File.ID)

# Output .rda file
usethis::use_data(psm_tmt_total)
