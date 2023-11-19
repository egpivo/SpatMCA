# SpatMCA Package

[![R build status](https://github.com/egpivo/SpatMCA/workflows/R-CMD-check/badge.svg)](https://github.com/egpivo/SpatMCA/actions)
[![Coverage Status](https://img.shields.io/codecov/c/github/egpivo/SpatMCA/master.svg)](https://codecov.io/github/egpivo/SpatMCA?branch=master)

### Description

**SpatMCA** is an R package designed for regularized maximum covariance analysis. It is a powerful tool for:

- Identifying smooth and localized couple patterns to understand how one spatial process affects another.
- Handling both regularly and irregularly spaced data, including 1D, 2D, and 3D datasets.
- Implementing the alternating direction method of multipliers (ADMM) algorithm.

### Installation

To install the current development version from GitHub:

```r
devtools::install_github("egpivo/SpatMCA")
```
Please note:
- Windows users require [Rtools](https://cran.r-project.org/bin/windows/Rtools/).

- Mac users need Xcode Command Line Tools and should install the library gfortran. You can do this by running the following commands in the terminal:
```bash
brew update
brew install gcc
```
For a detailed solution, refer to this [link](https://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/), or download and install the library [`gfortran`](https://github.com/fxcoudert/gfortran-for-macOS/releases) to resolve the "`ld: library not found for -lgfortran`" error.

### Usage
```r
library(SpatMCA)
spatmca(x1, x2, Y1, Y2, K = 1, num_cores = 1)
```
- Input: `x1`, `x2` (location matrices), `Y1`, `Y2` (data matrices), `K` (number of patterns), `num_cores` (number of CPU cores).
- Output: Provides information about the identified patterns.


### Author
 - [Wen-Ting Wang](https://www.linkedin.com/in/wen-ting-wang-6083a17b "Wen-Ting Wang")
 - [Hsin-Cheng Huang](https://sites.stat.sinica.edu.tw/hchuang/)
 
### Maintainer
[Wen-Ting Wang](https://www.linkedin.com/in/wen-ting-wang-6083a17b "Wen-Ting Wang")

### Reference
Wang, W.-T. and Huang, H.-C. (2017). [Regularized spatial maximum covariance analysis](https://arxiv.org/pdf/1705.02716.pdf), Environmetrics, 29, https://doi.org/10.1002/env.2481
 
### License
GPL-3
