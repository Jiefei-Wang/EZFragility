<!-- README.md is generated from README.Rmd. Please edit that file -->

# EZFragility: Epileptogenic Zone Localization Based on Fragility Index

[![](https://img.shields.io/badge/devel%20version-0.99.0-blue.svg)](https://github.com/Jiefei-Wang/EZFragility)
[![](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![](https://img.shields.io/github/languages/code-size/Jiefei-Wang/EZFragility.svg)](https://github.com/Jiefei-Wang/EZFragility)
[![](https://img.shields.io/github/last-commit/Jiefei-Wang/EZFragility.svg)](https://github.com/Jiefei-Wang/EZFragility/commits/main)
[![R-CMD-check](https://github.com/Jiefei-Wang/Fragility/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Jiefei-Wang/Fragility/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/Jiefei-Wang/Fragility/graph/badge.svg)](https://app.codecov.io/gh/Jiefei-Wang/Fragility)

To load the package

``` r
devtools::load_all()
```

The package contains an example data. To see it, type

``` r
pt01Epoch
```

For your test code, please consider creating a folder scripts and put
your code there. This folder will be ignored by git.

## TODO:

- Exported function names and parameters to snake case(e.g. nSearch -\>
  n_search)
- Unit test
- Vignette
- Check examples in the function documentation to make sure they are
  working
- Make sure all required functions/class have been exported
- Clear all error and warning in `devtools::check()` and
  `R CMD check --as-cran`
