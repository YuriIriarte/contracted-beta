[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20072138.svg)](https://doi.org/10.5281/zenodo.20072138)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

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
