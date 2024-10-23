## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
#install.packages("devtools")
#install.packages("knitr")
#install.packages("dglm")
#install.packages("doParallel")
#install.packages("foreach")
#install.packages("genio")
#install.packages("hglm")
#install.packages("Matrix")
#install.packages("onewaytests")
#library(devtools)
#devtools::install_github("kullrich/vGWAS", build_vignettes = TRUE, dependencies = FALSE)

## -----------------------------------------------------------------------------
library(vGWAS)

## -----------------------------------------------------------------------------
data(pheno)
data(geno)
data(chr)
data(map)

## -----------------------------------------------------------------------------
hist(pheno, breaks = 30, density = 15, col = "darkred")

## -----------------------------------------------------------------------------
table(chr)

## -----------------------------------------------------------------------------
data(geno)
dim(geno)
maf <- getMAF(geno.matrix = geno, geno.snp = "col")
geno.maf <- geno[, maf > 0.05]
dim(geno.maf)
chr.maf <- chr[maf > 0.05]
map.maf <- map[maf > 0.05]

## -----------------------------------------------------------------------------
vgwa <- vGWAS(
  phenotype = pheno,
  geno.matrix = geno,
  marker.map = map,
  chr.index = chr,
  pB = FALSE)

## -----------------------------------------------------------------------------
#vgwa <- vGWASparallel(
#  phenotype = pheno,
#  geno.matrix = geno,
#  geno.snp = "col",
#  marker.map = map,
#  chr.index = chr,
#  pB = FALSE,
#  ncores = 2)

## -----------------------------------------------------------------------------
plot(vgwa)

## -----------------------------------------------------------------------------
# get marker with lowest p.value
marker.lowest <- vgwa$p.value == min(vgwa$p.value)
vGWAS.variance(
  phenotype = pheno,
  marker.genotype = geno[, marker.lowest])

