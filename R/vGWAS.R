#' @title Variance Genome-wide Association
#' @name vGWAS
#' @aliases vGWAS
#' @description Variance Genome-wide association for using
#' nonparametric variance test
#' @usage vGWAS(phenotype, geno.matrix, kruskal.test = FALSE,
#' marker.map = NULL, chr.index = NULL, pB = TRUE)
#' @param phenotype a \code{numeric} or \code{logical} vector
#' of the phenotyic values. See \bold{Examples}.
#' @param geno.matrix a \code{matrix} or \code{data.frame} with
#' individuals as rows and markers as columns. The marker genotypes
#' for each marker are coded as one column. See \bold{Examples}.
#' @param kruskal.test a \code{logical} value specifying whether to use
#' Kruskal-Wallis statistic. The default option is \code{FALSE}, i.e.,
#' the usual ANOVA statistic is used in place of Kruskal-Wallis statistic.
#' @param marker.map a \code{numeric} vector giving the marker map
#' positions for each chromosome. See \bold{Examples}.
#' @param chr.index a \code{numeric} vector giving the chromosome index
#' for each marker. See \bold{Examples}.
#' @param pB show progress bar
#' @return a \code{data.frame} containing columns of \code{marker} names,
#' \code{chromosome} indices, \code{marker.map} positions,
#' test \code{statistic} values, and \code{p.value} for each position.
#' @seealso \code{\link{package-vGWAS}}
#' @references Shen, X., Pettersson, M., Ronnegard, L. and Carlborg, O.
#' (2011): \bold{Inheritance beyond plain heritability:
#' variance-controlling genes in \emph{Arabidopsis thaliana}}.
#' \emph{PLoS Genetics}, \bold{8}, e1002839.\cr
#' @references Ronnegard, L., Shen, X. and Alam, M. (2010):
#' \bold{hglm: A Package for Fitting Hierarchical Generalized
#' Linear Models}. \emph{The R Journal}, \bold{2}(2), 20-28.\cr
#' @examples
#' \donttest{
#' # ----- load data ----- #
#' data(pheno)
#' data(geno)
#' data(chr)
#' data(map)
#' # ----- variance GWA scan ----- #
#' vgwa <- vGWAS(phenotype = pheno, geno.matrix = geno,
#' marker.map = map, chr.index = chr, pb = FALSE)
#' # ----- visualize the scan ----- #
#' plot(vgwa)
#' summary(vgwa)
#' # ----- calculate the variance explained by the strongest marker ----- #
#' vGWAS.variance(phenotype = pheno,
#' marker.genotype = geno[,vgwa$p.value == min(vgwa$p.value)])
#' # ----- genomic control ----- #
#' vgwa2 <- vGWAS.gc(vgwa)
#' plot(vgwa2)
#' summary(vgwa2)
#' }
#' @importFrom stats anova lm median pchisq ppoints qchisq sd
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom graphics abline axis mtext plot points
#' @author Xia Shen
#' @export vGWAS
vGWAS <- function(
    phenotype,
    geno.matrix,
    kruskal.test = FALSE,
    marker.map = NULL,
    chr.index = NULL,
    pB = TRUE) {
    Call <- match.call()
    # ----- check phenotype ----- #
    if (!is.numeric(phenotype) & !is.logical(phenotype)) {
        stop('phenotype has to be numeric or logical.')
    }
    #if (heva)
    #{
    # cat('correcting phenotype using HEVA ...\n')
    # if (is.numeric(phenotype)) family <- gaussian() else family <- binomial()
    # phenotype <- vGWAS.heva(phenotype, geno.matrix, kinship, family = family)$corrected.phenotype
    # cat('phenotype corrected.\n')
    #}
    n <- length(phenotype)
    m <- ncol(geno.matrix)
    # ----- check genotypes ----- #
    if (!is.matrix(geno.matrix) & !is.data.frame(geno.matrix)) {
        stop('geno.matrix has to be a matrix or a data frame.')
    }
    # ----- check if data sizes match ----- #
    if (n != nrow(geno.matrix)) {
        stop('size of phenotype and geno.matrix do not match.')
    }
    if (!is.null(chr.index)) {
        if (m != length(chr.index)) {
            stop('size of chr.index and geno.matrix do not match.')
        }
    } else {
        chr.index <- rep(1, m)
    }
    if (!is.null(marker.map)) {
        if (m != length(marker.map)) {
            stop('size of marker.map and geno.matrix do not match.')
        }
    } else {
        tab.chr <- table(chr.index)
        marker.map <- c()
        for (i in 1:length(tab.chr)) {
          marker.map <- c(marker.map, 1:tab.chr[i])
        }
    }
    # ----- preallocation ----- #
    p.values <- statistics <- numeric(m)
    if (pB) {
        pb <- txtProgressBar(style = 3)
    }
    # ----- scan using Brown-Forsythe test -----#
    for (j in 1:m) {
        #  if (!heva) {
        test <- try(brown.forsythe.test(phenotype, as.factor(geno.matrix[,j]), kruskal.test = kruskal.test), silent  = TRUE)
        #  } else {
        #   test <- try(wilcox.test(phenotype ~ as.factor(geno.matrix[,j])), silent  = TRUE)
        #  }
        if (!inherits(test, 'try-error')) {
            p.values[j] <- test$p.value
            statistics[j] <- test$statistic
        } else {
            p.values[j] <- 1
            statistics[j] <- 0
        }
        if (pB) {
            setTxtProgressBar(pb, j/m)
        }
    }
    cat('\n')
    marker.names <- names(as.data.frame(geno.matrix))
    res <- data.frame(
        marker = marker.names,
        chromosome = chr.index,
        marker.map = marker.map,
        statistic = statistics,
        p.value = p.values)
    class(res) <- 'vGWAS'
    return(res)
}
