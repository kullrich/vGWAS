#' @title Brown-Forsythe's Test of Equality of Variances
#' @name bfmedian.test
#' @aliases brown.forsythe.test
#' @description The function performs the robust Brown-Forsythe
#' test using the group medians.
#' @usage bfmedian.test(y, group)
#' @details Levene (1960) proposed a test for homogeneity of variances
#' in \emph{k} groups which is based on the ANOVA statistic applied to
#' absolute deviations of observations from the corresponding group mean.
#' The robust Brown-Forsythe version of the Levene-type test substitutes
#' the group mean by the group median in the classical Levene statistic.
#' @param formula a formula of the form lhs ~ rhs where lhs gives the sample
#' values and rhs the corresponding groups.
#' @param data a tibble or data frame containing the variables in formula.
#' @param alpha the level of significance to assess the statistical difference.
#' Default is set to alpha = 0.05.
#' @param na.rm a logical value indicating whether NA values should be stripped
#' before the computation proceeds. Default us set to TRUE.
#' @param verbose a logical for printing output to R console.
#' @return A list with class "owt" containing the following components:
#' \item{statistic}{the Brown-Forsythe test statistic.}
#' \item{parameter}{the parameter(s) of the approximate F distribution of the
#' test statistic.}
#' \item{p.value}{the p-value of the test.}
#' \item{alpha}{the level of significance to assess the statistical difference.}
#' \item{method}{the character string "Brown-Forsythe-Median Test".}
#' \item{data}{a data frame containing the variables in which NA values
#' (if exist) are removed.}
#' \item{formula}{a formula of the form lhs ~ rhs where lhs gives the sample
#' values and rhs the corresponding groups.}
#' @note Modified from the \code{onewaytests} package and \code{vGWAS}.
#' @references Brown, M. B. and Forsythe, A.B. (1974). \bold{Robust tests for
#' equality of variances.}
#' \emph{Journal of the American Statistical Association}, \bold{69},
#' 364-367.\cr
#' @references Levene, H. (1960).
#' \bold{Robust Tests for Equality of Variances}, \emph{in Contributions
#' to Probability and Statistics}, ed. I. Olkin, Palo Alto, CA: Stanford Univ.
#' Press.\cr
#' @examples
#' data(pheno)
#' data(geno)
#' df <- data.frame(phenotype = pheno, genotype = as.factor(geno[, 911]))
#' bfmedian.test(phenotype ~ genotype, data = df)
#' @export bfmedian.test
#' @author Kristian Ullrich
#' @keywords htest
#' @importFrom stats anova lm median pchisq ppoints qchisq sd
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom graphics abline axis mtext plot points
bfmedian.test <- function (
    formula,
    data,
    alpha = 0.05,
    na.rm = TRUE,
    verbose = TRUE) {
    data <- model.frame(formula, data)
    dp <- as.character(formula)
    DNAME <- paste(dp[[2L]], "and", dp[[3L]])
    METHOD <- "Brown-Forsythe-Median Test based on deviations from the median"
    if (na.rm) {
        completeObs <- complete.cases(data)
        data <- data[completeObs, ]
    }
    if (any(colnames(data) == dp[[3L]]) == FALSE)
        stop("The name of group variable does not match the variable names in
             the data. The group variable must be one factor.")
    if (any(colnames(data) == dp[[2L]]) == FALSE)
        stop("The name of response variable does not match the variable names
             in the data.")
    y = data[[dp[[2L]]]]
    group = data[[dp[[3L]]]]
    if (!(is.factor(group) | is.character(group)))
        stop("The group variable must be a factor or a character.")
    if (is.character(group))
        group <- as.factor(group)
    if (!is.numeric(y))
        stop("The response must be a numeric variable.")
    n <- length(y)
    x.levels <- levels(factor(group))
    y.median = median(y)
    y.medians <- tapply(y, group, median)
    y.n <- tapply(y, group, length)
    # Calculate the absolute deviations from the medians ()
    y.deviations <- abs(y - y.medians[group])
    # Perform ANOVA on the absolute deviations
    lm_model <- lm(y.deviations ~ group)
    anova_result <- anova(lm_model)
    Ftest <- anova_result[1, "F value"]
    df1 <- anova_result[1, "Df"]
    df2 <- anova_result[2, "Df"]
    p.value <- anova_result[1, "Pr(>F)"]
    if (verbose) {
        cat("\n", "", METHOD, paste("(alpha = ", alpha, ")", sep = ""), "\n",
            sep = " ")
        cat("-------------------------------------------------------------",
            "\n", sep = " ")
        cat("  data :", DNAME, "\n\n", sep = " ")
        cat("  statistic  :", Ftest, "\n", sep = " ")
        cat("  num df     :", df1, "\n", sep = " ")
        cat("  denom df   :", df2, "\n", sep = " ")
        cat("  p.value    :", p.value, "\n\n", sep = " ")
        cat(
            if (p.value > alpha) {
                "  Result     : Difference is not statistically significant."
            }
            else {
                "  Result     : Difference is statistically significant."
            }, "\n")
        cat("-------------------------------------------------------------",
            "\n\n", sep = " ")
    }
    result <- list()
    result$statistic <- Ftest
    result$parameter <- c(df1, df2)
    result$p.value <- p.value
    result$alpha <- alpha
    result$method <- METHOD
    result$data <- data
    result$formula <- formula
    attr(result, "class") <- "owt"
    invisible(result)
}
