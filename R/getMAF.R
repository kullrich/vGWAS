#' @title Get minor-allele-frequency
#' @name getMAF
#' @aliases getMAF
#' @description Calculates minor-allele-frequency (MAF)
#' @usage getMAF(geno.matrix, geno.snp = "row", include.het = FALSE)
#' @param geno.matrix a \code{matrix} or \code{data.frame} or
#' \code{sparseMatrix} with individuals as columns and markers as rows
#' (geno.snp = "row") or individuals as rows and markers as columns
#' (geno.snp = "col").
#' @param geno.snp if individuals at columns and markers at rows use "row" else
#' if individuals at rows and markers at columns use "col"
#' @param include.het specify if heterozygous calls should be split and added
#' equally to homozygous ref and alt counts (default = FALSE)
#' @return a \code{vector} containing minor-allele-frequency values
#' @examples
#' data(geno)
#' maf <- getMAF(geno.matrix = geno, geno.snp = "col")
#' data(geno.num)
#' maf.num <- getMAF(geno.matrix = geno.num, geno.snp = "col")
#' data(geno.df)
#' maf.df <- getMAF(geno.matrix = geno.df, geno.snp = "row")
#' data(geno.sparse)
#' maf.sparse <- getMAF(geno.matrix = geno.sparse, geno.snp = "row")
#' @importFrom Matrix colSums rowSums
#' @author Kristian Ullrich
#' @export getMAF
getMAF<- function(
    geno.matrix,
    geno.snp = "row",
    include.het = FALSE) {
    # Internal functions
    getMAFnum <- function(snps_bin, include.het = FALSE) {
        # Tabulates values 0, 1, 2
        counts <- tabulate(as.numeric(snps_bin) + 1, nbins = 3)
        g0 <- counts[1]  # Homozygous 0 count
        g1 <- counts[2]  # Heterozygous 1 count
        g2 <- counts[3]  # Homozygous 2 count
        if (include.het) {
            g0 <- g0 + g1
            g2 <- g2 + g1
        }
        site_maf <- min(g0, g2) / (g0 + g2)
        return(site_maf)
    }
    getMAFsite <- function(snps_bin, include.het = FALSE) {
        counts <- table(snps_bin)
        site_maf <- min(counts) / sum(counts)
        return(site_maf)
    }
    # Check if geno.matrix is a matrix or a data frame
    if (!is.matrix(geno.matrix) &
        !is.data.frame(geno.matrix) &
        !is(geno.matrix, 'sparseMatrix')) {
        stop('geno.matrix has to be a matrix or a data frame or
             a sparse matrix.')
    }
    options(scipen = 22)
    if (is(geno.matrix, 'sparseMatrix')) {
        if (geno.snp == "row") {
            g0 <- Matrix::rowSums(geno.sparse == 0)  # Homozygous 0 count
            g1 <- Matrix::rowSums(geno.sparse == 1)  # Heterozygous 1 count
            g2 <- Matrix::rowSums(geno.sparse == 2)  # Homozygous 2 count
        } else {
            g0 <- Matrix::colums(geno.sparse == 0)  # Homozygous 0 count
            g1 <- Matrix::colSums(geno.sparse == 1)  # Heterozygous 1 count
            g2 <- Matrix::colSums(geno.sparse == 2)  # Homozygous 2 count
        }
        if (include.het) {
            g0 <- g0 + g1
            g2 <- g2 + g1
        }
        # Calculate MAF
        maf <- pmin(g0, g2) / (g0 + g2)
    } else if (is.matrix(geno.matrix)) {
        if (all(is.character(geno.matrix))) {
            if (geno.snp == "row") {
                maf <- apply(geno.matrix, 1, function(x) {
                    getMAFsite(x, include.het)})
            } else if (geno.snp == "col") {
                maf <- apply(geno.matrix, 2, function(x) {
                    getMAFsite(x, include.het)})
            }
        } else if (all(is.numeric(geno.matrix))) {
            if (geno.snp == "row") {
                maf <- apply(geno.matrix, 1, function(x) {
                    getMAFnum(x, include.het)})
            } else if (geno.snp == "col") {
                maf <- apply(geno.matrix, 2, function(x) {
                    getMAFnum(x, include.het)})
            }
        }
    } else {
        # geno.matrix is data.frame
        if (geno.snp == "row") {
            maf <- apply(geno.matrix, 1, function(x) {
                getMAFnum(x, include.het)})
        } else if (geno.snp == "col") {
            maf <- apply(geno.matrix, 2, function(x) {
                getMAFnum(x, include.het)})
        }
    }
    return(maf)
}
