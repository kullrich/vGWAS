#' @title Variance GWA Manhattan Plot
#' @name plot.vGWAS
#' @aliases plot.vGWAS
#' @description The function plots the variance GWA result for
#' the given scan object.
#' @usage \method{plot}{vGWAS}(x, sig.threshold = NULL, low.log.p = 0,
#' pch = 16, cex = 0.6, col.manhattan = c("slateblue4", "olivedrab"),
#' col.sig.threshold = "darkgoldenrod", ...)
#' @param x a result object from \code{vGWAS} scan.
#' It can be any \code{list} or \code{data.frame} that contains
#' \code{chromosome}, \code{marker.map}, and \code{p.value}, with
#' \code{class = 'vGWAS'}. See \code{\link{vGWAS}}.
#' @param sig.threshold a numeric value giving the significance
#' threshold for \code{-log(pvalues, 10)}. If \code{NULL},
#' Bonferroni correction will be used.
#' @param low.log.p a numeric value giving the lower limit of the
#' \code{-log(pvalues, 10)} to plot.
#' @param pch point character. See \code{\link{par}}.
#' @param cex size of points. See \code{\link{par}}.
#' @param col.manhattan two colors as a vector for the Manhattan plot.
#' @param col.sig.threshold one color for the significance threshold
#' @param ... not in use
#' @return a plot for viewing vGWAS result.
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
#' @export
plot.vGWAS <-
  function(x, sig.threshold = NULL, low.log.p = 0, pch = 16, cex = .6, col.manhattan = c('slateblue4', 'olivedrab'), col.sig.threshold = 'darkgoldenrod', ...)
  {
    tab.chr <- table(x$chromosome)
    chr <- as.numeric(names(tab.chr))
    ends <- cumsum(tab.chr)
    cumpos <- numeric(length(x$marker.map))
    cumpos[1:ends[1]] <- x$marker.map[1:ends[1]]
    logp <- -log(x$p.value, 10)
    if (length(chr) > 1)
    {
      for (i in 2:length(chr))
      {
        cumpos[(ends[i - 1] + 1):ends[i]] <- cumpos[ends[i - 1]] + x$marker.map[(ends[i - 1] + 1):ends[i]]
      }
    }
    if (is.null(sig.threshold))
    {
      sig.threshold <- -log(.05/length(x$marker.map), 10)
      cat('nominal significance threshold with Bonferroni correction for', length(x$marker.map), 'tests are calculated.\n')
    }
    cutp <- logp > low.log.p
    plot(cumpos, logp, type = 'n', ann = FALSE, axes = FALSE)
    at1 <- cumpos[1]
    points(cumpos[(1:ends[1])[cutp[1:ends[1]]]], logp[(1:ends[1])[cutp[1:ends[1]]]], pch = pch, cex = cex, col = col.manhattan[1])
    at1 <- c(at1, cumpos[ends[1]])
    at2 <- (cumpos[1] + cumpos[ends[1]])/2
    if (length(chr) > 1)
    {
      for (i in 2:length(chr))
      {
        points(cumpos[((ends[i - 1] + 1):ends[i])[cutp[(ends[i - 1] + 1):ends[i]]]], logp[((ends[i - 1] + 1):ends[i])[cutp[(ends[i - 1] + 1):ends[i]]]], pch = pch, cex = cex, col = col.manhattan[2 - i%%2])
        at1 <- c(at1, cumpos[ends[i]])
        at2 <- c(at2, (cumpos[ends[i - 1]] + cumpos[ends[i]])/2)
      }
    }
    abline(h = sig.threshold, lty = 2, col = col.sig.threshold)
    axis(1, at = at1, labels = FALSE)
    mtext(chr, 1, 0, at = at2)
    mtext('Chromosome', 1, 2)
    axis(2)
    mtext(expression(-log[10]~'('~italic(P)~-value~')'), 2, 3)
  }
