#' @title Variance GWA Summary
#' @name summary.vGWAS
#' @aliases summary.vGWAS
#' @description The function summarized the variance GWA
#' result for the given scan object.
#' @usage \method{summary}{vGWAS}(object, nrMarkers = 10, ...)
#' @param object a result object from \code{vGWAS} scan.
#' It can be any \code{list} or \code{data.frame} that contains
#' \code{chromosome}, \code{marker.map}, and \code{p.value},
#' with \code{class = 'vGWAS'}. See \code{\link{vGWAS}}.
#' @param nrMarkers a numeric value giving the number of top
#' markers to be summarized.
#' @param ... not in use
#' @return a summary for viewing vGWAS result.
#' @seealso \code{\link{package-vGWAS}}, \code{\link{vGWAS}}
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
#' @author Xia Shen
#' @importFrom stats anova lm median pchisq ppoints qchisq sd
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom graphics abline axis mtext plot points
#' @export
summary.vGWAS <- function(
    object,
    nrMarkers = 10,
    ...) {
    if(!class(object) == "vGWAS"){
        stop("data has to be of class: vGWAS")
    }
    pSort <- sort(object$p.value, index.return=T)
    topMarkers <- pSort$ix[1:nrMarkers]
    Pval <- object$p.value[topMarkers]
    chr <- object$chromosome[topMarkers]
    marker <- object$marker[topMarkers]
    map <- object$marker.map[topMarkers]
    result <- data.frame(marker, chr, map, Pval)
    print(paste("Top ", nrMarkers,
        " markers, sorted by p-value:", sep=""), quote=F)
    result
}
