#' @title Brown-Forsythe's Test of Equality of Variances
#' @name brown.forsythe.test
#' @aliases brown.forsythe.test
#' @description The function performs the robust Brown-Forsythe
#' test using the group medians. Instead of the ANOVA statistic,
#' the Kruskal-Wallis ANOVA may also be applied using this function.
#' @usage brown.forsythe.test(y, group, kruskal.test=FALSE)
#' @details Levene (1960) proposed a test for homogeneity of variances
#' in \emph{k} groups which is based on the ANOVA statistic applied to
#' absolute deviations of observations from the corresponding group mean.
#' The robust Brown-Forsythe version of the Levene-type test substites
#' the group mean by the group median in the classical Levene statistic.
#' @param y a numeric vector of data values.
#' @param group factor of the data.
#' @param kruskal.test a \code{logical} value specifying whether to use
#' Kruskal-Wallis statistic. The default option is \code{FALSE}, i.e.,
#' the usual ANOVA statistic is used in place of Kruskal-Wallis statistic.
#' @return A list with the following numeric components.
#' \item{statistic}{the value of the test statistic.}
#' \item{p.value}{the p-value of the test.}
#' \item{method}{type of test performed.}
#' \item{data.name}{a character string giving the name of the data.}
#' @note Modified from the \code{lawstat} package.
#' @references Brown, M. B. and Forsythe, A.B. (1974). \bold{Robust tests for equality of
#' variances.} \emph{Journal of the American Statistical Association}, \bold{69}, 364-367.\cr
#' @references Levene, H. (1960). \bold{Robust Tests for Equality of Variances}, \emph{in Contributions
#' to Probability and Statistics}, ed. I. Olkin, Palo Alto, CA: Stanford Univ. Press.\cr
#' @examples
#' data(pheno)
#' data(geno)
#' brown.forsythe.test(pheno, geno[,911])
#' @export brown.forsythe.test
#' @author Xia Shen
#' @keywords htest
#' @importFrom stats anova lm median pchisq ppoints qchisq sd
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom graphics abline axis mtext plot points
brown.forsythe.test <-
  function(y, group, kruskal.test = FALSE)
  {
    # ----- stop the code if the length of y does not match the length of group ----- #
    if (length(y) != length(group))
    {
      stop('The length of the data does not match the length of the group.')
    }
    # ----- assign stuffs ----- #
    DNAME <- deparse(substitute(y))
    y <- y[!is.na(y)]
    group <- group[!is.na(y)]
    # ----- sort the order just in case the input is not sorted by group ----- #
    reorder <- order(group)
    group <- group[reorder]
    y <- y[reorder]
    gr <- group
    group <- as.factor(group) # precautionary
    # ----- define the measure of central tendency (median) ----- #
    means <- tapply(y, group, median)
    METHOD <- 'Brown-Forsythe test based on the absolute deviations from the median'
    # ----- calculate the sample size of each group and absolute deviation from center ----- #
    n <- tapply(y, group, length)
    resp.mean <- abs(y - means[group])
    ngroup <- n[group]
    # ----- set d ----- #
    d <- group
    # ----- if the Kruskal-Wallis test is not used ----- #
    if (kruskal.test == FALSE)
    {
      statistic <- anova(lm(resp.mean ~ d))[1, 4]
      p.value <- anova(lm(resp.mean ~ d))[1, 5]
    }
    # ----- if the Kruskal-Wallis test is used ----- #
    else
    {
      METHOD <- paste('Rank-based (Kruskal-Wallis)', METHOD)
      ktest <- kruskal.test(resp.mean,d)
      statistic <- ktest$statistic
      p.value <- ktest$p.value
    }
    # ----- display output ----- #
    STATISTIC <- statistic
    names(STATISTIC) = 'Test Statistic'
    structure(list(statistic = STATISTIC, p.value = p.value, method = METHOD, data.name = DNAME), class = 'htest')
  }
