
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
This method implements an intracranial EEG (iEEG) marker of the
epileptogenic zone localization. In this method, seizures are
conceptualized as transitions from a stable networked system to an
unstable one. To quantify this, node fragility is computed from linear
network models, measuring each node’s susceptibility to destabilization.
There are significant details missing in (Li et al. 2017, 2021) to
reproduce the neural fragility method and adjust the parameters. This
Rpackage aims to identify and fill up the implementation details. It
will also allow users to test the method parameters on their data.

## EZFragility package tutorial

To load the package

``` r
devtools::load_all()
```

The package contains an example data. To see it, type

``` r
pt01Epochm1sp2s
```

The package contains an example results. To see it, type

``` r
pt01Fragm1sp2s
```

For your test code, please consider creating a folder scripts and put
your code there. This folder will be ignored by git. For explanations on
how to use the package please refer to the vignette
Intro_to_EZFragility.

## Implementation details

The method is based on building a discrete time linear system computing
a stable adjacency matrix A for the evolution of x(t).  
$x(t+1)=A x(t)$ with $x_i(t)$ the iEEG signal at time $t$ for electrode
$i$. A is computed for a series of time windows to derive the fragility
row.  
In this package, we are applying a ridge regression to fit the linear
operator A using the glmnet R library. In (Li et al. 2017, 2021), a
regularization parameter value of 1e-4 is recommended, however testing
on the data from patient pt01 from the Fragility data set (data subset
available in this package) this value does not ensure that A is always
stable. To tackle this issue, we have implemented a dichotomy to search
for the lowest stable lambda value rendering the matrix A stable (see R
function ridgesearchlambdadichomotomy in file ridge.r).

The method to compute the row perturbation is also not clear. To compute
the fragility row, a minimum 2-induced norm additive row perturbation
$\Delta$ is computed to destabilize the linear network placing an
eigenvalue of $A+\Delta$ at $\lambda=\sigma+j\omega$. The minimum norm
is a function of $\lambda$ given in (Li et al. 2017) (see function
fragilityRow in the scrip fragility.r), however the paper does not
describe how to choose $\lambda$ with $|\lambda|=1$. To tackle this
issue, we search for the value that minimize the norm of $\Delta$.

## TODO:

- Unit test (`test-main.R` from cecile)
- Vignette
- Make sure all required functions/class have been exported
- Clear all error and warning in `devtools::check()` and
  `R CMD check --as-cran`

Data:

- The data is too large. We need to find a way to reduce it.
  - Do we really need to recreate the result in the paper? For users,
    they just need to know how to use the package. If they want, they
    can download the data from the paper. Done

Ioannis:

- Branch1: Simplify Cecile’s `frag_stat` function in another branch (one
  or two branches depending on the amount of work) (Reviewer: Cecile)
- Branch2: Add hidden parameters to the Fragility class and (Reviewer:
  Jiefei)

Cecile:

- Branch3: Fix warning in `heatmap_frag` from ggplot2 and wording issue.
  (respect cammelCase `listelecmissing` -\> `listElecMissing`)
  (Reviewer: Jiefei)
- After Branch2: Use the hidden parameters in the Fragility class in
  visualization
- Make the epoch data within -1 to 2s
- why we call the data PT01? This number does not make any sense to
  users. Will we have PT02?
- What is example 3?

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
