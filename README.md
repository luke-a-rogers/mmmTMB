
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mmmTMB

[![Travis build
status](https://travis-ci.org/luke-a-rogers/mmmTMB.svg?branch=master)](https://travis-ci.org/luke-a-rogers/mmmTMB)
[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![Coverage
status](https://codecov.io/gh/luke-a-rogers/mmmTMB/branch/master/graph/badge.svg)](https://codecov.io/github/luke-a-rogers/mmmTMB?branch=master)

An R package to fit Markov movement models using TMB

## Installation

Installation and use of the `mmmTMB` package requires a [C++
compiler](https://support.rstudio.com/hc/en-us/articles/200486498-Package-Development-Prerequisites)
with `OpenMP` support.

``` r
devtools::install_github("luke-a-rogers/mmmTMB")
```

## Parallel processing

On my macOS machine running R 3.6.3 and following the R toolchain
install instructions from [The Coatless
Professor](https://thecoatlessprofessor.com/programming/cpp/r-compiler-tools-for-rcpp-on-macos-before-r-4.0.0/)
I made a new file `src/Makevars` with one line of code `PKG_CXXFLAGS =
$(SHLIB_OPENMP_CXXFLAGS)` followed by an empty line to make use of
parallel processing via OpenMP.
