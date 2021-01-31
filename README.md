[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/SpatMCA)](https://cran.r-project.org/package=SpatMCA)
[![Downloads](http://cranlogs.r-pkg.org/badges/grand-total/SpatMCA)](https://cran.r-project.org/package=SpatMCA)
[![Travis-CI Build Status](https://travis-ci.org/egpivo/SpatMCA.svg?branch=master)](https://travis-ci.org/egpivo/SpatMCA)
[![R build status](https://github.com/egpivo/SpatMCA/workflows/R-CMD-check/badge.svg)](https://github.com/egpivo/SpatMCA/actions)
[![Coverage Status](https://img.shields.io/codecov/c/github/egpivo/SpatMCA/master.svg)](https://codecov.io/github/egpivo/SpatMCA?branch=master)


# SpatMCA

### Description
***SpatMCA*** is an R package that provides regularized maximum covariance analysis, 

* identifying smooth and localized ***couple*** patterns to understand how one spatial process is affected by another
* suitable for either regularly or irregularly spaced data, including 1D, 2D, and 3D
* by the alternating direction method of multipliers (ADMM) algorithm

### Installation
To get the current released version from CRAN:

```r
install.packages("SpatMCA")
```

To get the current development version from GitHub:

```r
devtools::install_github("egpivo/SpatMCA")
```

To compile C++ code with the package [`RcppArmadillo`](https://cran.r-project.org/web/packages/RcppArmadillo/index.html),

 * Windows users require [Rtools](https://cran.r-project.org/bin/windows/Rtools/)
 * Mac users require Xcode Command Line Tools, and install the library gfortran by typing the following lines into terminal

  ```
   curl -O http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2
   sudo tar fvxz gfortran-4.8.2-darwin13.tar.bz2 -C /
  ```
  
More details can be found [here](http://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/).

### Author
 [Wen-Ting Wang](https://www.linkedin.com/in/wen-ting-wang-6083a17b "Wen-Ting Wang") and [Hsin-Cheng Huang](http://www.stat.sinica.edu.tw/hchuang/ "Hsin-Cheng Huang")
 
### Maintainer
[Wen-Ting Wang](https://www.linkedin.com/in/wen-ting-wang-6083a17b "Wen-Ting Wang")

### Reference
Wang, W.-T. and Huang, H.-C. (2017). [Regularized spatial maximum covaraince analysis](https://arxiv.org/pdf/1705.02716.pdf), Environmetrics, 29, https://doi.org/10.1002/env.2481
 
### License
  GPL-2
