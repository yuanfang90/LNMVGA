# LNMVGA
 This package provide a novel mixture of logistic-normal multinomial models for clustering microbiome data, or other compositional data. Parameter estimation are accomplished using variational Gaussian approximation to lift the computational overhead.

# How to Install

The preferred way to install this package is through through use devtools:

## Preparations

If you haven't installed the packages <span type="text/css">devtools</span> and <span type="text/css">roxygen2</span>, please do the following first:

```{r}
install.packages("devtools")
library("devtools")
devtools::install_github("klutometis/roxygen")
library(roxygen2)
```

## Install the package from GitHub
```r
devtools::install_github("yuanfang90/LNMVGA", upgrade_dependencies = FALSE)
```
