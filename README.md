
<!-- README.md is generated from README.Rmd. Please edit that file -->

# EZFragility: Epileptogenic Zone Localization Based on neural Fragility EEG marker

[![](https://img.shields.io/badge/devel%20version-0.99.0-blue.svg)](https://github.com/Jiefei-Wang/EZFragility)
[![](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![](https://img.shields.io/github/languages/code-size/Jiefei-Wang/EZFragility.svg)](https://github.com/Jiefei-Wang/EZFragility)
[![](https://img.shields.io/github/last-commit/Jiefei-Wang/EZFragility.svg)](https://github.com/Jiefei-Wang/EZFragility/commits/main)
[![R-CMD-check](https://github.com/Jiefei-Wang/Fragility/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Jiefei-Wang/Fragility/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/Jiefei-Wang/Fragility/graph/badge.svg)](https://app.codecov.io/gh/Jiefei-Wang/Fragility)

## Introduction

The goal of this Rpackage is to allow neuroscientists to reproduce and
test the neural fragility method described in (Li et al. 2017, 2021).
This method is a iEEG marker of the epileptogenic zone localization. In
this method, seizures are conceptualized as transitions from a stable
networked system to an unstable one. To quantify this, node fragility is
computed from linear network models, measuring each node’s
susceptibility to destabilization. There are significant details missing
in (Li et al. 2017, 2021) to reproduce the neural fragility method and
adjust the parameters. This Rpackage aims to identify and fill up the
implementation details. It will also allow users to test the method
parameters on their data.

## EZFragility package tutorial

To load the package

``` r
devtools::load_all()
```

The package contains an example data. To see it, type

``` r
pt01Epoch
```

The package contains an example results. To see it, type

``` r
pt01Frag
```

For your test code, please consider creating a folder scripts and put
your code there. This folder will be ignored by git.

## Implementation details

The method is based on building a discrete time linear system computing
a stable adjacency matrix A for the evolution of x(t). Equation … We
used a ridge regression … In … they recommend a lambda of 1e-4, however
testing on the data from pt01 (data available in this package) this
lambda value does not ensure that A is always stable.

The method to compute the row perturbation is also not clear

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

## References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-LiFragility2021" class="csl-entry">

Li, Adam, Chester Huynh, Zhary Fitzgerald, Iahn Cajigas, and Damina
Brusko. 2021. “Neural Fragility as an EEG Marker of the Seizure Onset
Zone.” *Nature Neuroscience* 24 (10): 1465–74.
<https://doi.org/10.1038/s41593-021-00901-w>.

</div>

<div id="ref-LiFragility2017" class="csl-entry">

Li, Adam, Sara Inati, Kareem Zaghloul, and Srivedi Sarma. 2017.
*Fragility in Epileptic Networks: The Epileptogenic Zone*. Lecture Notes
in Computer Science. IEEE. <https://doi.org/10.23919/ACC.2017.7963378>.

</div>

</div>
