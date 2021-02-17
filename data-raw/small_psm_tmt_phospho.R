# Read in PSM data
tmt_phospho <- read.delim(
  here::here("data-raw/Anja_phospho_lopit_phosphoproteome_rep2_PSMs.txt")
)

# Set seed for reproducible sampling
set.seed(2021)

# Randomly sample 5000 rows
tmt_phospho <- tmt_phospho[sample(nrow(tmt_phospho), 5000), ]

# Randomly sample 234 rows
tmt_phospho_filt <- tmt_phospho[sample(nrow(tmt_phospho), 234), ]

# Extract expression data (i.e. PSM abundances in each TMT channel)
abundance_cols <- colnames(tmt_phospho_filt)[grepl('Abundance.', colnames(tmt_phospho_filt))]
e_data <- as.matrix(tmt_phospho_filt[, abundance_cols])

# Extract feature data that we want (i.e. metadata about each PSM)
f_data <- tmt_phospho_filt[, c("Master.Protein.Accessions"), drop = FALSE]

# Construct MSnSet
small_psm_tmt_phospho <- MSnbase::MSnSet(exprs = e_data, fData = f_data)

# Output .rda file
save(small_psm_tmt_phospho, compress = "bzip2", file = here::here("inst/testdata/small_psm_tmt_phospho.rda"))
