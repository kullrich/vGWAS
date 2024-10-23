#' @title Variance Genome-wide Association (parallel)
#' @name vGWASparallel
#' @aliases vGWASparallel
#' @description Variance Genome-wide association for using
#' nonparametric variance test and other
#' @usage vGWASparallel(phenotype, geno.matrix, geno.snp = "row",
#' marker.map = NULL, chr.index = NULL, method = "bfmedian",
#' p.adjust.method = "none", include.het = FALSE, pB = TRUE, ncores = 1)
#' @param phenotype a \code{numeric} or \code{logical} vector
#' of the phenotyic values.
#' @param geno.matrix a \code{matrix} or \code{data.frame} or
#' \code{sparseMatrix} with individuals as columns and markers as rows
#' (geno.snp = "row") or individuals as rows and markers as columns
#' (geno.snp = "col").
#' @param marker.map a \code{numeric} vector giving the marker map
#' positions for each chromosome.
#' @param chr.index a \code{numeric} vector giving the chromosome index
#' for each marker.
#' @param geno.snp if individuals at columns and markers at rows use "row" else
#' if individuals at rows and markers at columns use "col"
#' @param method the test method to use (default = bfmedian).
#' Default is set to the Brown-Forsythe's Test of Equality of Variances
#' using group medians.
#' There are 31 other tests available via the onewaytests package:
#' Alvandi's F test ("af"),
#' Alexander-Govern test ("ag"),
#' Alvandi's generalized p-value ("agp"),
#' One-way analysis of variance ("aov"),
#' Approximate F test ("ap"),
#' Adjusted Welch's heteroscedastic F test ("aw"),
#' B square test ("b2"),
#' Brown-Forsythe test ("bf"),
#' Box F test ("box"),
#' Cochran test ("cochran"),
#' Generalized tests equivalent to Parametric Bootstrap ("gtb"),
#' Generalized tests equivalent to Fiducial tests ("gtf"),
#' Variance homogeneity tests ("homog"),
#' James second order test ("james"),
#' Johansen F test ("johansen"),
#' Kruskal-Wallis test ("kw"),
#' Modified Brown-Forsythe test ("mbf"),
#' Mann-Whitney U test ("mw"),
#' Anderson-Darling normility test ("nor_ad"),
#' Cramer-vin Mises normility test ("nor_cvm"),
#' Kolmogorov-Smirnov normility test ("nor_ks"),
#' Pearson Chi-square normility test ("nor_pct"),
#' Shapiro-Wilk normility test ("nor_sw"),
#' Shapiro-Francia normility test ("nor_sf"),
#' Permutation F test ("pf"),
#' Scott-Smith test ("ss"),
#' Student's t-test ("st"),
#' Welch-Aspin test ("wa"),
#' Welch's heteroscedastic F test with trimmed means and Winsorized variances
#' ("welch"),
#' Weerahandi's generalized F test ("wgf"),
#' Welch's t-test ("wt").
#' @param test.alpa the level of significance to assess the statistical
#' difference. Default is set to alpha = 0.05.
#' @param test.na.rm a logical value indicating whether NA values should be
#' stripped before the computation proceeds. Default us set to TRUE.
#' @param p.adjust.method correction method (default = "none").
#' There are 8 p-value correction methods available via the p.adjust function:
#' "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
#' @param include.het specify if heterozygous calls should be split and added
#' equally to homozygous ref and alt counts (default = FALSE)
#' @param pB show progress bar
#' @param ncores number of cores to parallelize (default = 1)
#' @return a \code{data.frame} containing columns of \code{marker} names,
#' \code{chromosome} indices, \code{marker.map} positions,
#' test \code{statistic} values, and \code{p.value} for each position.
#' @seealso \code{\link{package-vGWAS}} \code{onewaytests}
#' @references Shen, X., Pettersson, M., Ronnegard, L. and Carlborg, O.
#' (2011): \bold{Inheritance beyond plain heritability:
#' variance-controlling genes in \emph{Arabidopsis thaliana}}.
#' \emph{PLoS Genetics}, \bold{8}, e1002839.\cr
#' @references Ronnegard, L., Shen, X. and Alam, M. (2010):
#' \bold{hglm: A Package for Fitting Hierarchical Generalized
#' Linear Models}. \emph{The R Journal}, \bold{2}(2), 20-28.\cr
#' @examples
#' # ----- load data ----- #
#' data(pheno)
#' data(geno)
#' data(chr)
#' data(map)
#' # ----- variance GWA scan ----- #
#' vgwa <- vGWASparallel(phenotype = pheno, geno.matrix = geno,
#' marker.map = map, chr.index = chr,
#' geno.snp = "col", pb = FALSE)
#' # ----- other test GWA scan ----- #
#' vgwa.mw <- vGWASparallel(phenotype = pheno, geno.matrix = geno,
#' marker.map = map, chr.index = chr,
#' geno.snp = "col", method = "mw", pb = FALSE)
#' # ----- multiple cores ----- #
#' vgwa.st <- vGWASparallel(phenotype = pheno, geno.matrix = geno,
#' marker.map = map, chr.index = chr,
#' geno.snp = "col", method = "st", ncores = 2, pb = FALSE)
#' @importFrom stats anova lm median pchisq ppoints qchisq sd p.adjust
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom graphics abline axis mtext plot points
#' @importFrom Matrix Matrix
#' @importFrom onewaytests onewaytests
#' @importFrom parallel stopCluster makeCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @author Xia Shen
#' @author Kristian Ullrich
#' @export vGWASparallel
vGWASparallel <- function(
    phenotype,
    geno.matrix,
    marker.map = NULL,
    chr.index = NULL,
    geno.snp = "row",
    method = "bfmedian",
    test.alpha = 0.05,
    test.na.rm = TRUE,
    p.adjust.method = "none",
    include.het = FALSE,
    pB = TRUE,
    ncores = 1) {
    Call <- match.call()
    # ----- check phenotype ----- #
    if (!is.numeric(phenotype) & !is.logical(phenotype)) {
        stop('phenotype has to be numeric or logical.')
    }
    n <- length(phenotype)
    # ----- check genotypes ----- #
    # Check if geno.matrix is a matrix or a data frame
    if (!is.matrix(geno.matrix) &
        !is.data.frame(geno.matrix) &
        !is(geno.matrix, 'sparseMatrix')) {
      stop('geno.matrix has to be a matrix or a data frame or a sparse matrix.')
    }
    if (geno.snp == "row") {
        m <- nrow(geno.matrix)
    } else if (geno.snp == "col") {
        m <- ncol(geno.matrix)
    }
    # ----- check if data sizes match ----- #
    # Ensure the dimensions match between phenotype and geno.matrix
    if (geno.snp == "row") {
        if (n != ncol(geno.matrix)) {
            stop('Size of phenotype and geno.matrix do not match.')
        }
    } else if (geno.snp == "col") {
        if (n != nrow(geno.matrix)) {
            stop('Size of phenotype and geno.matrix do not match.')
        }
    }
    # Check chr.index dimension
    if (!is.null(chr.index)) {
        if (m != length(chr.index)) {
            stop('Size of chr.index and geno.matrix do not match.')
        }
    } else {
        chr.index <- rep(1, m)
    }
    # Check marker.map dimension
    if (!is.null(marker.map)) {
        if (m != length(marker.map)) {
            stop('Size of marker.map and geno.matrix do not match.')
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
    # ----- scan -----#
    if (ncores > 1) {
        if(.Platform$OS.type == "windows"){
            cl <- parallel::makeCluster(ncores)
        } else {
            cl <- parallel::makeForkCluster(ncores)
        }
        chunk_indices <- split(seq(from = 1, to = m),
                               cut(seq(from = 1, to = m), ncores))
        doParallel::registerDoParallel(cl)
        results <- foreach::foreach(batch = chunk_indices,
                                    .combine = 'rbind') %dopar% {
            local_p_values <- numeric(length(batch))
            local_statistics <- numeric(length(batch))
            if (geno.snp == "row") {
                geno.matrix.batch <- geno.matrix[batch, ]
                bm <- nrow(geno.matrix.batch)
            } else if (geno.snp == "col") {
                geno.matrix.batch <- geno.matrix[, batch]
                bm <- ncol(geno.matrix.batch)
            }
            for (j in 1:bm) {
            # Process by row or column based on geno.snp
            if (geno.snp == "row") {
                data <- data.frame(phenotype = phenotype,
                                   genotype = as.factor(geno.matrix.batch[j, ]))
            } else if (geno.snp == "col") {
                data <- data.frame(phenotype = phenotype,
                                   genotype = as.factor(geno.matrix.batch[, j]))
            }
            if (method == "af") {
                test <- try(onewaytests::af.test(phenotype ~ genotype,
                    data = data,
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "ag") {
                test <- try(onewaytests::ag.test(phenotype ~ genotype,
                    data = data,
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "agp") {
                test <- try(onewaytests::agp.test(phenotype ~ genotype,
                    data = data,
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "aov") {
                test <- try(onewaytests::aov.test(phenotype ~ genotype,
                    data = data,
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "ap") {
                test <- try(onewaytests::ap.test(phenotype ~ genotype,
                    data = data,
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "aw") {
                test <- try(onewaytests::aw.test(phenotype ~ genotype,
                    data = data,
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "b2") {
                test <- try(onewaytests::b2.test(phenotype ~ genotype,
                    data = data,
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "bf") {
                test <- try(onewaytests::bf.test(phenotype ~ genotype,
                    data = data,
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "bfmedian") {
                test <- try(bfmedian.test(phenotype ~ genotype,
                    data = data,
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "box") {
                test <- try(onewaytests::box.test(phenotype ~ genotype,
                    data = data,
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "cochran") {
                test <- try(onewaytests::cochran.test(phenotype ~ genotype,
                    data = data,
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "gtb") {
                test <- try(onewaytests::gp.test(phenotype ~ genotype,
                    data = data,
                    method = "gtb",
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "gtf") {
                test <- try(onewaytests::gp.test(phenotype ~ genotype,
                    data = data,
                    method = "gtf",
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "homog") {
                test <- try(onewaytests::homog.test(phenotype ~ genotype,
                    data = data,
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "james") {
                test <- try(onewaytests::james.test(phenotype ~ genotype,
                    data = data,
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "johansen") {
                test <- try(onewaytests::johansen.test(phenotype ~ genotype,
                    data = data,
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "kw") {
                test <- try(onewaytests::kw.test(phenotype ~ genotype,
                    data = data,
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "mbf") {
                test <- try(onewaytests::mbf.test(phenotype ~ genotype,
                    data = data,
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "mw") {
                test <- try(onewaytests::mw.test(phenotype ~ genotype,
                    data = data,
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "nor_ad") {
                test <- try(onewaytests::nor.test(phenotype ~ genotype,
                    data = data,
                    method = "AD",
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "nor_cvm") {
                test <- try(onewaytests::nor.test(phenotype ~ genotype,
                    data = data,
                    method = "CVM",
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "nor_ks") {
                test <- try(onewaytests::nor.test(phenotype ~ genotype,
                    data = data,
                    method = "LT",
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "nor_pct") {
                test <- try(onewaytests::nor.test(phenotype ~ genotype,
                    data = data,
                    method = "PT",
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "nor_sw") {
                test <- try(onewaytests::nor.test(phenotype ~ genotype,
                    data = data,
                    method = "SW",
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "nor_sf") {
                test <- try(onewaytests::nor.test(phenotype ~ genotype,
                    data = data,
                    method = "SF",
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "pf") {
                test <- try(onewaytests::pf.test(phenotype ~ genotype,
                    data = data,
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "ss") {
                test <- try(onewaytests::ss.test(phenotype ~ genotype,
                    data = data,
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "st") {
                test <- try(onewaytests::st.test(phenotype ~ genotype,
                    data = data,
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "wa") {
                test <- try(onewaytests::wa.test(phenotype ~ genotype,
                    data = data,
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "welch") {
                test <- try(onewaytests::welch.test(phenotype ~ genotype,
                    data = data,
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "wgf") {
                test <- try(onewaytests::wgf.test(phenotype ~ genotype,
                    data = data,
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "wt") {
                test <- try(onewaytests::wt.test(phenotype ~ genotype,
                    data = data,
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else {
                stop('Wrong test method.')
            }
            # Store the results
            if (!inherits(test, 'try-error')) {
                local_p_values[j] <- ifelse(is.na(test$p.value), 1,
                                            test$p.value)
                local_statistics[j] <- ifelse(is.na(test$statistic), 0,
                                              test$statistic)
            } else {
                local_p_values[j] <- 1
                local_statistics[j] <- 0
            }
        }
        # Combine the local results into a list and return
        return(data.frame(p.values = local_p_values,
                          statistics = local_statistics, indices = batch))
      }
      # Collect results from all batches
      p.values <- unlist(results[["p.values"]][order(results[["indices"]])])
      statistics <- unlist(results[["statistics"]][order(results[["indices"]])])
      parallel::stopCluster(cl)
    } else {
        for (j in 1:m) {
            # Create a data frame with the phenotype and the genotype at SNP j
            if (geno.snp == "row") {
                data <- data.frame(phenotype = phenotype,
                                   genotype = as.factor(geno.matrix[j, ]))
            } else if (geno.snp == "col") {
                data <- data.frame(phenotype = phenotype,
                                   genotype = as.factor(geno.matrix[, j]))
            }
            if (method == "af") {
                test <- try(onewaytests::af.test(phenotype ~ genotype,
                    data = data,
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "ag") {
                test <- try(onewaytests::ag.test(phenotype ~ genotype,
                    data = data,
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "agp") {
                test <- try(onewaytests::agp.test(phenotype ~ genotype,
                    data = data,
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "aov") {
                test <- try(onewaytests::aov.test(phenotype ~ genotype,
                    data = data,
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "ap") {
                test <- try(onewaytests::ap.test(phenotype ~ genotype,
                    data = data,
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "aw") {
                test <- try(onewaytests::aw.test(phenotype ~ genotype,
                    data = data,
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "b2") {
                test <- try(onewaytests::b2.test(phenotype ~ genotype,
                    data = data,
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "bf") {
                test <- try(onewaytests::bf.test(phenotype ~ genotype,
                    data = data,
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "bfmedian") {
                test <- try(bfmedian.test(phenotype ~ genotype,
                    data = data,
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "box") {
                test <- try(onewaytests::box.test(phenotype ~ genotype,
                    data = data,
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "cochran") {
                test <- try(onewaytests::cochran.test(phenotype ~ genotype,
                    data = data,
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "gtb") {
                test <- try(onewaytests::gp.test(phenotype ~ genotype,
                    data = data,
                    method = "gtb",
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "gtf") {
                test <- try(onewaytests::gp.test(phenotype ~ genotype,
                    data = data,
                    method = "gtf",
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "homog") {
                test <- try(onewaytests::homog.test(phenotype ~ genotype,
                    data = data,
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "james") {
                test <- try(onewaytests::james.test(phenotype ~ genotype,
                    data = data,
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "johansen") {
                test <- try(onewaytests::johansen.test(phenotype ~ genotype,
                    data = data,
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "kw") {
                test <- try(onewaytests::kw.test(phenotype ~ genotype,
                    data = data,
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "mbf") {
                test <- try(onewaytests::mbf.test(phenotype ~ genotype,
                    data = data,
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "mw") {
                test <- try(onewaytests::mw.test(phenotype ~ genotype,
                    data = data,
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "nor_ad") {
                test <- try(onewaytests::nor.test(phenotype ~ genotype,
                    data = data,
                    method = "AD",
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "nor_cvm") {
                test <- try(onewaytests::nor.test(phenotype ~ genotype,
                    data = data,
                    method = "CVM",
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "nor_ks") {
                test <- try(onewaytests::nor.test(phenotype ~ genotype,
                    data = data,
                    method = "LT",
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "nor_pct") {
                test <- try(onewaytests::nor.test(phenotype ~ genotype,
                    data = data,
                    method = "PT",
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "nor_sw") {
                test <- try(onewaytests::nor.test(phenotype ~ genotype,
                    data = data,
                    method = "SW",
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "nor_sf") {
                test <- try(onewaytests::nor.test(phenotype ~ genotype,
                    data = data,
                    method = "SF",
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "pf") {
                test <- try(onewaytests::pf.test(phenotype ~ genotype,
                    data = data,
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "ss") {
                test <- try(onewaytests::ss.test(phenotype ~ genotype,
                    data = data,
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "st") {
                test <- try(onewaytests::st.test(phenotype ~ genotype,
                    data = data,
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "wa") {
                test <- try(onewaytests::wa.test(phenotype ~ genotype,
                    data = data,
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "welch") {
                test <- try(onewaytests::welch.test(phenotype ~ genotype,
                    data = data,
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "wgf") {
                test <- try(onewaytests::wgf.test(phenotype ~ genotype,
                    data = data,
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else if (method == "wt") {
                test <- try(onewaytests::wt.test(phenotype ~ genotype,
                    data = data,
                    alpha = test.alpha,
                    na.rm = test.na.rm,
                    verbose = FALSE), silent = TRUE)
            } else {
                stop('Wrong test method.')
            }
            if (!inherits(test, 'try-error')) {
                p.values[j] <- ifelse(is.na(test$p.value), 1, test$p.value)
                statistics[j] <- ifelse(is.na(test$statistic), 0,
                                        test$statistic)
            } else {
                p.values[j] <- 1
                statistics[j] <- 0
            }
            if (pB) {
                setTxtProgressBar(pb, j / m)
            }
        }
    }
    if (pB) {
        cat('\n')
    }
    if (geno.snp == "row") {
        if (!is(geno.matrix, 'sparseMatrix')) {
            marker.names <- rownames(as.data.frame(geno.matrix))
        } else {
            marker.names <- rownames(geno.matrix)
        }
    } else if (geno.snp == "col") {
        if (!is(geno.matrix, 'sparseMatrix')) {
            marker.names <- colnames(as.data.frame(geno.matrix))
        } else {
            marker.names <- colnames(geno.matrix)
        }
    }
    res <- data.frame(
        marker = marker.names,
        chromosome = chr.index,
        marker.map = marker.map,
        statistic = statistics,
        p.value = p.values
    )
    if (p.adjust.method != "none") {
        res[["p.value"]] <- stats::p.adjust(res[["p.value"]],
            method = p.adjust.method)
    }
    class(res) <- 'vGWAS'
    return(res)
}
