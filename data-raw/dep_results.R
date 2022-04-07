library(DEP)
library(dplyr)

# Read in data from DEP package
proteins <- DEP::UbiLength
expdesign <- DEP::UbiLength_ExpDesign

# Perform differential enrichment analysis
results <- LFQ(proteins, expdesign, fun = "MinProb", type = "control",
               control = "Ctrl", alpha = 0.05, lfc = 1)$results

# Split the UniProt accessions in results table into 100 accessions chunks
# (Otherwise the UniProt API will give a 400 error)
accessions_chunks <- split(results$ID, ceiling(seq_along(results$ID) / 100))

# Map UniProt accessions to Ensembl genome IDs
ensembl_lst <- vector(mode = "list", length = length(accessions_chunks))

for (i in seq_along(accessions_chunks)) {
  payload <- list(
    query = paste(accessions_chunks[[i]], collapse = " "),
    from = "ACC+ID",
    to = "ENSEMBL_ID",
    format = "tab"
  )

  response <- httr::GET(
    url = "https://www.uniprot.org/uploadlists/",
    query = payload
  )
  message(response$status_code)

  ensembl_lst[[i]] <- read.table(text = httr::content(response),
                                 sep = "\t", header = TRUE,
                                 col.names = c("ID", "ensembl"))
}

ensembl_mappings <- Reduce(rbind, ensembl_lst) %>%
  distinct(ID, .keep_all = TRUE) %>%
  distinct(ensembl, .keep_all = TRUE)

# Get lengths for each protein in results table
protein_lst <- vector(mode = "list", length = length(accessions_chunks))

for (i in seq_along(accessions_chunks)) {
  payload <- list(
    query = paste(accessions_chunks[[i]], collapse = " "),
    from = "ACC+ID",
    to = "ACC",
    format = "tab",
    columns = paste(c("id", "length"), collapse = ",")
  )

  response <- httr::GET(
    url = "https://www.uniprot.org/uploadlists/",
    query = payload
  )
  message(response$status_code)

  protein_lst[[i]] <- read.table(text = httr::content(response),
                                 sep = "\t", header = TRUE,
                                 col.names = c("ID_short", "length", "ID"))

}

protein_lengths <- Reduce(rbind, protein_lst) %>%
  distinct(ID, .keep_all = TRUE)

# Join 3 data.frames together
# (keep only rows with valid length and Ensembl gene ID)
# (ensure no duplicates in: ensembl, ID, and ID_short)
dep_results <- left_join(ensembl_mappings, results, by = "ID") %>%
  left_join(protein_lengths, by = "ID") %>%
  select(name, ID, ID_short, ensembl, length, everything()) %>%
  distinct(ID_short, .keep_all = TRUE)

# Output .rda file
usethis::use_data(dep_results)
