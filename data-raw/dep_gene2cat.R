library(camprotR)
library(dplyr)
# also requires goseq package

# keep only rows with unique ensembl IDs
dep_filt <- dep_results %>%
  distinct(ensembl, .keep_all = TRUE)

# obtain GO terms annotated to these proteins (just for test purposes)
dep_gene2cat <- goseq::getgo(dep_filt$ensembl, genome = "hg19", id = "ensGene",
                             fetch.cats=c("GO:CC", "GO:BP", "GO:MF"))

# Output .rda file
usethis::use_data(dep_gene2cat, compress = TRUE)
