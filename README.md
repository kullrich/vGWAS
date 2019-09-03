# vGWAS
[![Build Status](https://travis-ci.com/kullrich/vGWAS.svg?branch=master)](https://travis-ci.com/kullrich/vGWAS)

Variance Heterogeneity Genome-wide Association Study - Reimplementation
=========
This repository is a reimplementation from the original `vGWAS` R package from [Xia Shen](https://github.com/xiashen).

see the original publication
[Inheritance Beyond Plain Heritability: Variance-Controlling Genes in Arabidopsis thaliana](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1002839)

## Installation

### R specific installation prerequisites

```
install.packages("devtools")
install.packages("knitr")
install.packages("hglm")
install.packages("dglm")
```

Third install `vGWAS` package from [github](https://github.com/kullrich) using the [devtools](https://cran.r-project.org/web/packages/devtools/index.html) package.

```
library(devtools)
install_github("kullrich/vGWAS", build_vignettes = TRUE, dependencies = FALSE)
```

### Vignettes

These vignettes introduce `distIUPAC`

- [01. vGWAS Tutorial](https://github.com/kullrich/vGWAS/blob/master/vignettes/vGWAS.Rmd)

```
library(vGWAS)
vignette("vGWAS")
