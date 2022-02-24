library(dplyr)
# files from https://github.com/CambridgeCentreForProteomics/infSILAC/tree/master/raw/exp_2
silac_pep_infiles <- Sys.glob(here::here('data-raw/OOPS*_PeptideGroups.txt'))
silac_psm_infiles <- Sys.glob(here::here('data-raw/OOPS*_PSMs.txt'))

pep_data <- lapply(silac_pep_infiles, read.delim)
psm_data <- lapply(silac_psm_infiles, read.delim)

# Sample 100 proteins to keep

set.seed(0)

proteins <- lapply(pep_data, function(x) x$Master.Protein.Accessions) %>%
  Reduce(f='intersect') %>%
  unique()

keep_proteins <- sample(proteins[!grepl(';', proteins)], 100)

filter_inf <- function(x){
  x %>%
    filter((grepl('cRAP', Master.Protein.Accessions) |
              Master.Protein.Accessions %in% keep_proteins))
}

psm_oops_1 <- psm_data[[1]] %>% filter_inf()
psm_oops_2 <- psm_data[[2]] %>% filter_inf()
psm_oops_3 <- psm_data[[3]] %>% filter_inf()
psm_oops_4 <- psm_data[[4]] %>% filter_inf()
pep_oops_1 <- pep_data[[1]] %>% filter_inf()
pep_oops_2 <- pep_data[[2]] %>% filter_inf()
pep_oops_3 <- pep_data[[3]] %>% filter_inf()
pep_oops_4 <- pep_data[[4]] %>% filter_inf()

usethis::use_data(psm_oops_1, overwrite=TRUE)
usethis::use_data(psm_oops_2, overwrite=TRUE)
usethis::use_data(psm_oops_3, overwrite=TRUE)
usethis::use_data(psm_oops_4, overwrite=TRUE)
usethis::use_data(pep_oops_1, overwrite=TRUE)
usethis::use_data(pep_oops_2, overwrite=TRUE)
usethis::use_data(pep_oops_3, overwrite=TRUE)
usethis::use_data(pep_oops_4, overwrite=TRUE)

