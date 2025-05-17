# vGWAS
[![Build Status](https://travis-ci.com/kullrich/vGWAS.svg?branch=master)](https://travis-ci.com/kullrich/vGWAS)

Variance Heterogeneity Genome-wide Association Study - Reimplementation
=========

R package source code: https://github.com/kullrich/vGWAS

R package pages: https://kullrich.github.io/vGWAS/

R package issues: https://github.com/kullrich/vGWAS/issues

This repository is a reimplementation from the original `vGWAS` R package from [Xia Shen](https://github.com/xiashen).

see the original publication
[Inheritance Beyond Plain Heritability: Variance-Controlling Genes in Arabidopsis thaliana](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1002839)

some function has been added to directly perform GWAS on genotype data obtained
via [plink](https://www.cog-genomics.org/plink/)

the function `vGWASparallel` has been added to perform statistical tests in parallel
and work on genotype data encoded as sparse matrix

## Installation

### R specific installation prerequisites

```
install.packages("devtools")
install.packages("knitr")
install.packages("dglm")
install.packages("doParallel")
install.packages("foreach")
install.packages("genio")
install.packages("hglm")
install.packages("Matrix")
install.packages("onewaytests")
```

Install `vGWAS` package from [github](https://github.com/kullrich) using the [devtools](https://cran.r-project.org/web/packages/devtools/index.html) package.

```
library(devtools)
devtools::install_github("kullrich/vGWAS", build_vignettes = TRUE, dependencies = FALSE)
```

## Quick start

```
library(vGWAS)
data(pheno)
data(geno.sparse)
data(chr)
data(map)
vgwa <- vGWASparallel(
  phenotype = pheno,
  geno.matrix = geno.sparse,
  marker.map = map,
  chr.index = chr,
  geno.snp = "row"
)
plot(vgwa)
```

## Vignettes

These vignettes introduce `vGWAS`

- [01. vGWAS basic tutorial](https://github.com/kullrich/vGWAS/blob/master/vignettes/vGWAS.Rmd)
- [02. vGWAS plink tutorial](https://github.com/kullrich/vGWAS/blob/master/vignettes/vGWASsparse.Rmd)


## Code of Conduct - Participation guidelines

This repository adhere to [Contributor Covenant](http://contributor-covenant.org) code of conduct for in any interactions you have within this project. (see [Code of Conduct](https://github.com/kullrich/CRBHits/-/blob/devel/CODE_OF_CONDUCT.md))

See also the policy against sexualized discrimination, harassment and violence for the Max Planck Society [Code-of-Conduct](https://www.mpg.de/11961177/code-of-conduct-en.pdf).

By contributing to this project, you agree to abide by its terms.
