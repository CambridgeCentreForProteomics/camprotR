camprotR
================

<!-- badges: start -->

[![R-CMD-check](https://github.com/CambridgeCentreForProteomics/camprotR/workflows/check-bioc/badge.svg)](https://github.com/CambridgeCentreForProteomics/camprotR/actions)
[![codecov](https://codecov.io/gh/CambridgeCentreForProteomics/camprotR/branch/master/graph/badge.svg)](https://codecov.io/gh/CambridgeCentreForProteomics/camprotR?branch=master)
[![License:
MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

The `camprotR` package is an internal package written and used by the
[Cambridge Centre for Proteomics](https://proteomics.bio.cam.ac.uk) to
facilitate our proteomic analyses. It is opinionated by design and most
functions assume the data has been processed by Proteome Discoverer.

## Install

``` r
remotes::install_github("CambridgeCentreForProteomics/camprotR", build_vignettes = TRUE, dependencies='Suggests')
```

## How to use

Several HTML vignettes are installed with this package which go through
how to use different aspects of the package. You can see a list of
available vignettes with:

``` r
browseVignettes("camprotR")
```
