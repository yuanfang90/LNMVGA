# LNMVGA
 This package provide a novel mixture of logistic-normal multinomial models for clustering microbiome data, or other compositional data. Parameter estimation are accomplished using variational Gaussian approximation to lift the computational overhead.

## How to Install

The preferred way to install this package is through using `devtools::install_github()`:

### Preparations

If you haven't installed the package [devtools](https://cran.r-project.org/web/packages/devtools/index.html), please do the following first:

```{r}
install.packages("devtools")
library("devtools")
```

### Install the package from GitHub

```r
devtools::install_github("yuanfang90/LNMVGA", upgrade_dependencies = FALSE)
```
### Then...enjoy using the package

```{r}
library(LNMVGA)
```
