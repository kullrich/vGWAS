#' @title Genomic Control for vGWAS
#' @name vGWAS.gc
#' @aliases vGWAS.gc
#' @description The function does genomic control for the
#' variance GWA result object.
#' @usage vGWAS.gc(object, plot = TRUE, proportion = 1, ...)
#' @param object a result object from \code{vGWAS} scan.
#' It can be any \code{list} or \code{data.frame} that contains
#' \code{chromosome}, \code{marker.map}, and \code{p.value},
#' with \code{class = 'vGWAS'}. See \code{\link{vGWAS}}.
#' @param plot a logical value turning on/off the QQ plot
#' for genomic control.
#' @param proportion a numeric value between 0 and 1 giving
#' the proportion of obtained p-values to be used for genomic control.
#' @param ... not in use
#' @return
#' \item{lambda}{estimated inflation ratio.}
#' \item{lambda.se}{standard error of the estimated inflation ratio.}
#' \item{gc.p.value}{p-values after genomic control.}
#' @seealso \code{\link{package-vGWAS}}, \code{\link{vGWAS}}
#' @references Shen, X., Pettersson, M., Ronnegard, L. and Carlborg, O.
#' (2011): \bold{Inheritance beyond plain heritability:
#' variance-controlling genes in \emph{Arabidopsis thaliana}}.
#' \emph{PLoS Genetics}, \bold{8}, e1002839.\cr
#' @examples
#' \donttest{
#' # ----- load data ----- #
#' data(pheno)
#' data(geno)
#' data(chr)
#' data(map)
#' # ----- variance GWA scan ----- #
#' vgwa <- vGWAS(phenotype = pheno, geno.matrix = geno,
#' marker.map = map, chr.index = chr)
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
#' @export vGWAS.gc
vGWAS.gc <- function(object, plot = TRUE, proportion = 1, ...)
{
  if (proportion > 1 || proportion <= 0)
    stop('proportion argument should be greater then zero and less than or equal to one.')
  ntp <- round(proportion * length(object$p.value))
  if (ntp <= 1)
    stop('too few valid measurments.')
  if (ntp < 10)
    warning(paste('number of points is fairly small:', ntp))
  if (min(object$p.value) < 0)
    stop('data argument has values <0')
  if (max(object$p.value) <= 1) {
    data <- data0 <- qchisq(object$p.value, 1, lower.tail = FALSE)
  }
  data <- sort(data)
  ppoi <- ppoints(data)
  ppoi <- sort(qchisq(1 - ppoi, 1))
  data <- data[1:ntp]
  ppoi <- ppoi[1:ntp]
  s <- summary(lm(data ~ 0 + ppoi))$coeff
  out <- object
  out$lambda <- s[1, 1]
  out$lambda.se <- s[1, 2]
  if (plot) {
    lim <- c(0, max(data, ppoi, na.rm = TRUE))
    plot(ppoi, data, xlab = expression('Expected'~chi^2), ylab = expression('Observed'~chi^2), ...)
    abline(a = 0, b = 1, col = 4, lwd = 2.4)
    if (out$lambda > 1) co <- 2 else co <- 3
    abline(a = 0, b = (s[1, 1]), col = co, lwd = 2.4)
  }
  out$p.value <- pchisq(data0/out$lambda, 1, lower.tail = FALSE)
  return(out)
}
