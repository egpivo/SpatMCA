# SpatMCA: Regularized Spatial Maximum Covariance Analysis

[![License](https://eddelbuettel.github.io/badges/GPL2+.svg)](https://www.gnu.org/licenses/gpl-2.0.html)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/SpatMCA?color=green)](https://cran.r-project.org/package=SpatMCA)
[![R build status](https://github.com/egpivo/SpatMCA/workflows/R-CMD-check/badge.svg)](https://github.com/egpivo/SpatMCA/actions)
[![Coverage Status](https://img.shields.io/codecov/c/github/egpivo/SpatMCA/master.svg)](https://codecov.io/github/egpivo/SpatMCA?branch=master)
[![Downloads (monthly)](https://cranlogs.r-pkg.org/badges/SpatMCA?color=brightgreen)](https://www.r-pkg.org/pkg/SpatMCA)
[![Downloads (total)](https://cranlogs.r-pkg.org/badges/grand-total/SpatMCA?color=brightgreen)](https://www.r-pkg.org/pkg/SpatMCA)
[![Environmetrics](https://img.shields.io/badge/Environmetrics-10.1002%2Fenv.2481-brightgreen)](https://doi.org/10.1002/env.2481)

## Description

**SpatMCA** is an R package designed for regularized maximum covariance analysis. It serves as a powerful tool for:

- Identifying smooth and localized coupling patterns to understand how one spatial process affects another.
- Handling both regularly and irregularly spaced data, spanning 1D, 2D, and 3D datasets.
- Implementing the alternating direction method of multipliers (ADMM) algorithm.


## Installation
You can install the **SpatMCA** package using one of the following methods:

### Install from CRAN:
```r
install.packages("SpatMCA")
```

### Install the current development version from GitHub:
```r
remotes::install_github("egpivo/SpatMCA")
```
#### Please Note:
- **Windows Users:** Ensure that you have [Rtools](https://cran.r-project.org/bin/windows/Rtools/) installed before proceeding with the installation.

- **Mac Users:** You need Xcode Command Line Tools and should install the library [`gfortran`](https://github.com/fxcoudert/gfortran-for-macOS/releases). Follow these steps in the terminal:
    ```bash
    brew update
    brew install gcc
    ```
    For a detailed solution, refer to this [link](https://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/), or download and install the library [`gfortran`](https://github.com/fxcoudert/gfortran-for-macOS/releases) to resolve the "`ld: library not found for -lgfortran`" error.


### Usage
To perform regularized maximum covariance analysis using **SpatMCA**, follow these steps:

```r
library(SpatMCA)
spatmca(x1, x2, Y1, Y2, K = 1, num_cores = 1)
```
#### Parameters:
  - `x1`, `x2`: Location matrices.
  - `Y1`, `Y2`: Data matrices.
  - `K`: Number of patterns.
  - `num_cores`: Number of CPU cores.
#### Output:
Provides information about the identified patterns

## Authors
 - [Wen-Ting Wang](https://www.linkedin.com/in/wtwang) ([GitHub](https://www.github.com/egpivo))
 - [Hsin-Cheng Huang](https://sites.stat.sinica.edu.tw/hchuang/)
 
## Maintainer
[Wen-Ting Wang](https://www.linkedin.com/in/wtwang) ([GitHub](https://www.github.com/egpivo))

## Reference
Wang, W.-T. and Huang, H.-C. (2018). [Regularized spatial maximum covariance analysis](https://arxiv.org/pdf/1705.02716.pdf), Environmetrics, 29, https://doi.org/10.1002/env.2481
 
## License
GPL (>= 2)

## Citation
1. To cite package ‘SpatMCA’ in publications use:
```
  Wang W, Huang H (2023). _SpatMCA: Regularized Spatial Maximum Covariance Analysis_.
  R package version 1.0.2.6, <https://github.com/egpivo/SpatMCA>.
```
2. A BibTeX entry for LaTeX users is
```
  @Manual{,
    title = {SpatMCA: Regularized Spatial Maximum Covariance Analysis},
    author = {Wen-Ting Wang and Hsin-Cheng Huang},
    year = {2023},
    note = {R package version 1.0.2.6},
    url = {https://github.com/egpivo/SpatMCA},
  }
```
