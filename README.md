
<!-- README.md is generated from README.Rmd. Please edit that file -->

# selEnergyPermR

<!-- badges: start -->
<!-- badges: end -->

Hinton, A.L. & Mucha, P. J., A simultaneous feature selection and
compositional association test for detecting sparse associations in
high-dimensional metagenomic data. 12 July 2021, PREPRINT (Version 1)
available at Research Square
<https://doi.org/10.21203/rs.3.rs-703177/v1>

The goal of selEnergyPermR package is to provide an easy to use set of
reusable functions for performing the selEnergyPerm methods for sparse
association testing in compositional data.

## Prerequisites

In order to use the ‘selEnergyPermR’ package the ‘diffCompVarRcpp’
package must first be installed. This package contains the cpp
implementation of the differential compositonal variation scoring
algorithm required for efficiently ranking logratios.

``` r
install.packages("devtools")
devtools::install_github(repo = "andrew84830813/diffCompVarRcpp",
                         dependencies = T)

```

## Installation

Next install the latest version of selEnergyPermR from github with:

``` r
devtools::install_github(repo = "andrew84830813/selEnergyPermR",
                         dependencies = T)
```

## Example

This is a basic example which demonstrate how to use ‘selEnergyPermR’
for testing associations in count compositional data:

``` r
library(diffCompVarRcpp)
library(selEnergyPermR)

## basic example code
```

Load the example dataset and view the first 5 rows and columns

``` r
data("Uguanda_PIH")

(Uguanda_PIH[1:5,1:5])
#>   Status
#> 1    pos
#> 2    neg
#> 3    pos
#> 4    pos
#> 5    pos
#>   Abditibacteriota.Abditibacteria.Abditibacteriales.Abditibacteriaceae.Abditibacterium
#> 1                                                                                    0
#> 2                                                                                    0
#> 3                                                                                    0
#> 4                                                                                    0
#> 5                                                                                    0
#>   Actinobacteriota.Acidimicrobiia.null.null.null
#> 1                                              0
#> 2                                              0
#> 3                                              0
#> 4                                              0
#> 5                                              0
#>   Actinobacteriota.Actinobacteria.Actinomycetales.Actinomycetaceae.Actinomyces
#> 1                                                                            0
#> 2                                                                            0
#> 3                                                                            0
#> 4                                                                            0
#> 5                                                                            0
#>   Actinobacteriota.Actinobacteria.Actinomycetales.Actinomycetaceae.F0332
#> 1                                                                      0
#> 2                                                                      0
#> 3                                                                      0
#> 4                                                                      0
#> 5                                                                      0
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/master/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
