# vGWAS
[![Build Status](https://travis-ci.com/kullrich/vGWAS.svg?branch=master)](https://travis-ci.com/kullrich/vGWAS)

Variance Heterogeneity Genome-wide Association Study - Reimplementation
=========
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
