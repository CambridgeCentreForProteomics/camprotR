% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/count_features_per_protein.R
\name{count_features_per_protein}
\alias{count_features_per_protein}
\title{Count the number of features (rows) per protein}
\usage{
count_features_per_protein(obj, master_prot_col = "Master.Protein.Accessions")
}
\arguments{
\item{obj}{\code{MSnSet}. Contains PSMs or peptides.}

\item{master_prot_col}{\code{character} Column name for master protein ID}
}
\value{
\code{tibble} with feature counts per sample per protein.
}
\description{
For differential testing with
\url{https://www.bioconductor.org/packages/release/bioc/html/DEqMS.html}{DEqMS},
one needs to identify the number of features, e.g PSMs or peptides per protein.
This function returns the number of features per protein per sample.
}
\examples{
# Use a small example TMT dataset included with the camprotR package
df <- psm_tmt_total

# Make an MSnSet
df_exprs <- as.matrix(df[, grep("Abundance", colnames(df))])
colnames(df_exprs) <- gsub("Abundance\\\\.", "", colnames(df_exprs))

df_fData <- df[, grep("Abundance", colnames(df), invert = TRUE)]

psm <- MSnbase::MSnSet(exprs = df_exprs, fData = df_fData)

# Count the number of PSMs per protein
count_features_per_protein(psm, master_prot_col = "Master.Protein.Accessions")

}
