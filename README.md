[![R build status](https://github.com/egpivo/SpatMCA/workflows/R-CMD-check/badge.svg)](https://github.com/egpivo/SpatMCA/actions)
[![Coverage Status](https://img.shields.io/codecov/c/github/egpivo/SpatMCA/master.svg)](https://codecov.io/github/egpivo/SpatMCA?branch=master)


# SpatMCA

### Description
***SpatMCA*** is an R package that provides regularized maximum covariance analysis, 

* identifying smooth and localized ***couple*** patterns to understand how one spatial process is affected by another
* suitable for either regularly or irregularly spaced data, including 1D, 2D, and 3D
* by the alternating direction method of multipliers (ADMM) algorithm

### Installation
To get the current development version from GitHub:

```r
devtools::install_github("egpivo/SpatMCA")
```

To compile C++ code with the package [`RcppArmadillo`](https://cran.r-project.org/web/packages/RcppArmadillo/index.html),

 * Windows users require [Rtools](https://cran.r-project.org/bin/windows/Rtools/)
 * Mac users require Xcode Command Line Tools, and install the library gfortran by typing the following lines into terminal
    ```
    brew update
    brew install gcc
    ```
    The detailed solution is describd [here](https://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/), or download and install the library [`gfortran`](https://github.com/fxcoudert/gfortran-for-macOS/releases) to solve the error `ld: library not found for -lgfortran`.

### Author
 [Wen-Ting Wang](https://www.linkedin.com/in/wen-ting-wang-6083a17b "Wen-Ting Wang") and [Hsin-Cheng Huang](https://sites.stat.sinica.edu.tw/hchuang/)
 
### Maintainer
[Wen-Ting Wang](https://www.linkedin.com/in/wen-ting-wang-6083a17b "Wen-Ting Wang")

### Reference
Wang, W.-T. and Huang, H.-C. (2017). [Regularized spatial maximum covariance analysis](https://arxiv.org/pdf/1705.02716.pdf), Environmetrics, 29, https://doi.org/10.1002/env.2481
 
### License
  GPL-2
