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

