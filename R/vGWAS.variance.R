#' @title Calculating Variance Explained by A Single Marker
#' @name vGWAS.variance
#' @aliases vGWAS.variance
#' @description The function calculates and reports the variance
#' explained for a single marker by fitting a double generalized
#' linear model. It gives both the variance explained by the mean
#' and variance parts of model.
#' @usage vGWAS.variance(phenotype, marker.genotype, print = TRUE)
#' @param phenotype a \code{numeric} vector of the phenotyic values.
#' See \bold{Examples}.
#' @param marker.genotype a \code{numeric} or \code{character} or
#' \code{factor} vector of the genotypes of a single marker.
#' See \bold{Examples}.
#' @param print a \code{logical} value. If \code{FALSE},
#' the heritability values will be returned for storage.
#' @details The \bold{Value} will only be available if
#' \code{only.print = FALSE}.
#' @return
#' \item{variance.mean}{the variance explained by the mean part of model.}
#' \item{variance.disp}{the variance explained by the variance part of model.}
#' @seealso \code{\link{package-vGWAS}}, \code{\link{vGWAS}}, \code{\link{plot.vGWAS}}
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
#' }
#' @author Xia Shen
#' @importFrom stats anova lm median pchisq ppoints qchisq sd
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom graphics abline axis mtext plot points
#' @export vGWAS.variance
vGWAS.variance <-
  function(phenotype, marker.genotype, print = TRUE)
  {
    # ----- check phenotype ----- #
    if (!is.numeric(phenotype))
    {
      stop('phenotype has to be numeric.')
    }
    # ----- check marker genotype ----- #
    if (!is.numeric(marker.genotype) & !is.character(marker.genotype) & !is.factor(marker.genotype))
    {
      stop('marker genotype has a wrong format.')
    }
    # ----- check if data sizes match ----- #
    if (length(phenotype) != length(marker.genotype))
    {
      stop('size of phenotype and marker genotype do not match.')
    }
    #dm <- dglm(phenotype ~ as.factor(marker.genotype), ~ as.factor(marker.genotype))
    #df.fd <- length(levels(as.factor(marker.genotype))) - 1
    #p.mean = summary(dm)$coef[2,4]
    #fd <- '~ 1'
    #fd <- parse(text = fd)
    #mode(fd) <- 'call'
    #lik <- rep(dm$m2loglik, 2)
    #names(lik) <- c('Mean', 'Full')
    #ncall <- dm$call
    #ncall['dformula'] <- fd
    #lik['Mean'] <- eval(ncall)$m2loglik
    #LRT = as.numeric(lik['Mean'] - lik['Full'])
    #p.disp = pchisq(LRT, df.fd, lower.tail = FALSE) ###
    #heritability.mean <- (dm$null.deviance - dm$deviance)/dm$null.deviance
    #heritability.disp0 <- (dm$dispersion.fit$null.deviance - dm$dispersion.fit$deviance)/dm$dispersion.fit$null.deviance
    #heritability.disp <- heritability.disp0*(1 - heritability.mean)
    tab <- table(marker.genotype)
    genos <- names(tab)
    if (length(genos) != 2) stop('Incorrect number of genotypes for calculating variance explained.')
    y1 <- phenotype[marker.genotype == genos[1]]
    y2 <- phenotype[marker.genotype == genos[2]]
    mu1 <- mean(y1)
    mu2 <- mean(y2)
    s1 <- sd(y1)
    s2 <- sd(y2)
    p <- tab[1]/length(phenotype)
    vp <- p*s1**2 + (1 - p)*s2**2 + p*(1 - p)*(mu1 - mu2)**2
    vm <- p*(1 - p)*(mu1 - mu2)**2
    vv <- p*(1 - p)*(s1 - s2)**2
    ve <- (p*s1 + (1 - p)*s2)**2
    if (print) {
      cat('variance explained by the mean part of model:\n')
      cat(round(vm/vp*100, digits = 2), '%\n')
      #cat(round(vm*100, digits = 2), '%, p-value =', p.mean, '\n')
      cat('variance explained by the variance part of model:\n')
      cat(round(vv/vp*100, digits = 2), '%\n')
      #cat(round(vv*100, digits = 2), '%, p-value =', p.disp, '\n')
      cat('variance explained in total:\n')
      cat(round((vm + vv)/vp*100, digits = 2), '%\n')
    }
    return(list(vm = vm, vv = vv, ve = ve, vp = vp))
    #p.mean = p.mean, p.variance = p.disp, LRT.statistic.variance = LRT, df.variance = df.fd
  }
