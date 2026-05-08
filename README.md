# CB: An R package for modeling unit data using the Contracted-Beta distribution (v1.0.1)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20072138.svg)](https://doi.org/10.5281/zenodo.20072138)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

R package for modeling data on the unit interval $(0,1)$ using the Contracted-Beta (CB) distribution.

## Overview

The **CB** package implements the Contracted-Beta distribution, a flexible model for bounded data on the unit interval $(0,1)$.

The package provides tools for:

- Density, distribution, quantile, and random generation functions  
- Moment-based and maximum likelihood estimation  
- Descriptive measures (mean, variance, skewness, kurtosis)  
- Model comparison and goodness-of-fit procedures  

All scripts used in the associated manuscript (simulations, figures, and real-data applications) are included for reproducibility.

## Installation

You can install the stable version used in the paper directly from GitHub:

```r
install.packages("remotes")
remotes::install_github("YuriIriarte/contracted-beta")
```

## Basic usage

```markdown
## Basic Usage

```r
library(CB)

# Generate data
x <- rCB(100, alpha = 2, beta = 3)

# Fit model (MLE)
fit <- fitCB_mle(x)

# Summary
fit$par
fit$AIC
```

## Reproducibility Resources

All scripts used in the simulation studies, figures, and real-data
applications reported in the paper *Flexible Modeling for Unit Data
via Random Contraction Mechanisms: Inference, Simulation, and
Applications* are distributed within the package.

After installing and loading the package, the available scripts can
be listed using:

```r
list.files(system.file("paper", package = "CB"))
```

For example, the Monte Carlo simulation script can be accessed
directly through:

```r
file.edit(system.file("paper/simulations.R", package = "CB"))
```

Additional scripts for figures and applications can be accessed in
the same way.

Documentation for all implemented functions is available through the
standard R help system. For example:

```r
?dCB
?fitCB_mle
?gofCB_boot
```

## Citation

If you use this package, please cite:

```r
@Manual{,
  title = {CB: R package for modeling unit data with the Contracted-Beta distribution},
  author = {Yuri A. Iriarte},
  year = {2026},
  note = {R package version 1.0.0},
  doi = {10.5281/zenodo.20072138},
  url = {https://doi.org/10.5281/zenodo.20072138},
}
```
Iriarte, Y. A. (2026). CB: R package for modeling unit data using the Contracted-Beta distribution.
DOI: https://doi.org/10.5281/zenodo.20072138
