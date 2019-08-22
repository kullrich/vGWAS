## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
#install.packages("devtools")
#install.packages("hglm")
#install.packages("dglm")
#library(devtools)
#install_github("kullrich/vGWAS", build_vignettes = TRUE, dependencies = FALSE)

## ------------------------------------------------------------------------
library(vGWAS)
#vignette("vGWAS")

## ------------------------------------------------------------------------
data(pheno)
data(geno)
data(chr)
data(map)

## ------------------------------------------------------------------------
hist(pheno, breaks = 30, density = 15, col = "darkred")

## ------------------------------------------------------------------------
table(chr)

## ------------------------------------------------------------------------
vgwa <- vGWAS(phenotype = pheno, geno.matrix = geno, marker.map = map, chr.index = chr, pB = FALSE)
#vgwa <- vGWAS(phenotype = pheno, geno.matrix = geno, marker.map = map, chr.index = chr, pB = TRUE)

## ------------------------------------------------------------------------
plot(vgwa)

## ------------------------------------------------------------------------
vGWAS.variance(phenotype = pheno, marker.genotype = geno[,vgwa$p.value == min(vgwa$p.value)])

## ------------------------------------------------------------------------
geno.coding <- matrix(0, nrow(geno), ncol(geno))
#pb <- txtProgressBar(style = 3)
for (j in 1:ncol(geno)) {
       geno.coding[,j] <- as.numeric(geno[,j] == names(table(geno[,j]))[1])*2 - 1
       #setTxtProgressBar(pb, j/ncol(geno))
}
image(tcrossprod(geno.coding))

## ------------------------------------------------------------------------
#vgwa2 <- vGWAS(phenotype = pheno, geno.matrix = geno.coding, heva = TRUE, marker.map = map, chr.index = chr)
vgwa2 <- vGWAS.gc(vgwa)

## ------------------------------------------------------------------------
plot(vgwa2)

