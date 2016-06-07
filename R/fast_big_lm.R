## bigLm.R: Rcpp/Eigen/bigmemory implementation of lm()
##
## Copyright (C)  2011 - 2015  Douglas Bates, Dirk Eddelbuettel and Romain Francois
##
## Copyright (C)  2016  Jared Huling
##
## This file is based off of code from RcppEigen.
##
## RcppEigen is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.
##
## RcppEigen is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with RcppEigen.  If not, see <http://www.gnu.org/licenses/>.



#' fast and memory efficient linear model fitting
#'
#' @param X input model matrix. must be a big.matrix object (type = 8 for double)
#' @param y numeric response vector of length nobs.
#' @param method an integer scalar with value 0 for the LLT Cholesky or 1 for the LDLT Cholesky
#' @return A list with the elements
#' \item{coefficients}{a vector of coefficients}
#' \item{se}{a vector of the standard errors of the coefficient estimates}
#' \item{rank}{a scalar denoting the computed rank of the model matrix}
#' \item{df.residual}{a scalar denoting the degrees of freedom in the model}
#' \item{residuals}{the vector of residuals}
#' \item{s}{a numeric scalar - the root mean square for residuals}
#' \item{fitted.values}{the vector of fitted values}
#' @export
#' @examples
#'
#' library(bigmemory)
#'
#' nrows <- 50000
#' ncols <- 50
#' bkFile <- "bigmat2.bk"
#' descFile <- "bigmatk2.desc"
#' bigmat <- filebacked.big.matrix(nrow=nrows, ncol=ncols, type="double",
#'                                 backingfile=bkFile, backingpath=".",
#'                                 descriptorfile=descFile,
#'                                 dimnames=c(NULL,NULL))
#'
#' # Each column value with be the column number multiplied by
#' # samples from a standard normal distribution.
#' set.seed(123)
#' for (i in 1:ncols) bigmat[,i] = rnorm(nrows)*i
#'
#' y <- rnorm(nrows) + bigmat[,1]
#'
#' system.time(lmr1 <- bigLmPure(bigmat, y))
#'
#' system.time(lmr2 <- lm.fit(x = bigmat[,], y = y))
#'
#' max(abs(coef(lmr1) - coef(lmr2)))
#'
#'
bigLmPure <- function(X, y, method = 0L) {

    stopifnot( is.big.matrix(X), is.numeric(y), NROW(y) == nrow(X) )

    .Call("bigLm_Impl", X@address, y, method, colnames(X), PACKAGE="bigFastlm")
}

#' fast and memory efficient linear model fitting
#'
#' @param X input model matrix. must be a big.matrix object (type = 8 for double)
#' @param y numeric response vector of length nobs.
#' @param method an integer scalar with value 0 for the LLT Cholesky or 1 for the LDLT Cholesky
#' @return A list object with S3 class "bigLm" with the elements
#' \item{coefficients}{a vector of coefficients}
#' \item{se}{a vector of the standard errors of the coefficient estimates}
#' \item{rank}{a scalar denoting the computed rank of the model matrix}
#' \item{df.residual}{a scalar denoting the degrees of freedom in the model}
#' \item{residuals}{the vector of residuals}
#' \item{s}{a numeric scalar - the root mean square for residuals}
#' \item{fitted.values}{the vector of fitted values}
#' @useDynLib bigFastlm
#' @import Rcpp
#' @export
#' @examples
#'
#' library(bigmemory)
#'
#' nrows <- 50000
#' ncols <- 50
#' bkFile <- "bigmat.bk"
#' descFile <- "bigmatk.desc"
#' bigmat <- filebacked.big.matrix(nrow=nrows, ncol=ncols, type="double",
#'                                 backingfile=bkFile, backingpath=".",
#'                                 descriptorfile=descFile,
#'                                 dimnames=c(NULL,NULL))
#'
#' # Each column value with be the column number multiplied by
#' # samples from a standard normal distribution.
#' set.seed(123)
#' for (i in 1:ncols) bigmat[,i] = rnorm(nrows)*i
#'
#' y <- rnorm(nrows) + bigmat[,1]
#'
#' system.time(lmr1 <- bigLm(bigmat, y))
#'
#' system.time(lmr2 <- lm.fit(x = bigmat[,], y = y))
#'
#' max(abs(coef(lmr1) - coef(lmr2)))
#'
#'
bigLm <- function(X, ...) UseMethod("bigLm")



#' summary method for bigLm fitted objects
#'
#' @param object bigLm fitted object
#' @param ... not used
#' @return a summary.bigLm object
#' @rdname summary
#' @method summary bigLm
#' @export
#' @examples
#'
#' library(bigmemory)
#'
#' nrows <- 50000
#' ncols <- 15
#' bkFile <- "bigmat4.bk"
#' descFile <- "bigmatk4.desc"
#' bigmat <- filebacked.big.matrix(nrow=nrows, ncol=ncols, type="double",
#'                                 backingfile=bkFile, backingpath=".",
#'                                 descriptorfile=descFile,
#'                                 dimnames=c(NULL,NULL))
#'
#' # Each column value with be the column number multiplied by
#' # samples from a standard normal distribution.
#' set.seed(123)
#' for (i in 1:ncols) bigmat[,i] = rnorm(nrows)*i
#'
#' y <- rnorm(nrows) + bigmat[,1] - bigmat[,2]
#'
#' system.time(lmr1 <- bigLm(bigmat, y))
#'
#' summary(lmr1)
#'
#'
summary.bigLm <- function(object, ...)
{
    coef <- object$coefficients
    se   <- object$se
    tval <- coef/se

    object$coefficients <- cbind(Estimate     = coef,
                                 "Std. Error" = se,
                                 "t value"    = tval,
                                 "Pr(>|t|)"   = 2*pt(-abs(tval), df=object$df))

    ## cf src/stats/R/lm.R and case with no weights and an intercept
    f <- object$fitted.values
    r <- object$residuals
    #mss <- sum((f - mean(f))^2)
    mss <- if (object$intercept) sum((f - mean(f))^2) else sum(f^2)
    rss <- sum(r^2)

    object$r.squared <- mss/(mss + rss)
    df.int <- if (object$intercept) 1L else 0L
    n <- length(f)
    rdf <- object$df
    object$adj.r.squared <- 1 - (1 - object$r.squared) * ((n - df.int)/rdf)
    class(object) <- "summary.bigLm"
    object
}


#' bigLm default
#'
#' @param ... not used
#' @rdname bigLm
#' @method bigLm default
#' @export
bigLm.default <- function(X, y, method = 0L, ...) {

    y <- as.numeric(y)

    res <- bigLmPure(X, y, method)
    res$call <- match.call()
    res$intercept <- any(big.colMax(X) == big.colMin(X))

    class(res) <- "bigLm"
    res
}

#' print method for bigLm objects
#'
#' @rdname print
#' @method print bigLm
#' @export
print.bigLm <- function(x, ...) {
    cat("\nCall:\n")
    print(x$call)
    cat("\nCoefficients:\n")
    print(x$coefficients, digits=5)
}

#' print method for summary.bigLm objects
#'
#' @param x a "summary.bigLm" object
#' @param ... not used
#' @rdname print
#' @method print summary.bigLm
#' @export
print.summary.bigLm <- function(x, ...) {
    cat("\nCall:\n")
    print(x$call)
    cat("\nResiduals:\n")
    digits <- max(3, getOption("digits") - 3)
    print(summary(x$residuals, digits=digits)[-4])
    cat("\n")

    printCoefmat(x$coefficients, P.values=TRUE, has.Pvalue=TRUE, ...)
    cat("\nResidual standard error: ", formatC(x$s, digits=digits), " on ",
        formatC(x$df), " degrees of freedom\n", sep="")
    cat("Multiple R-squared: ", formatC(x$r.squared, digits=digits),
        ",\tAdjusted R-squared: ",formatC(x$adj.r.squared, digits=digits),
        "\n", sep="")
    invisible(x)
}

#' Prediction method for bigLm fitted objects
#'
#' @param object fitted "bigLm" model object
#' @param newdata big.matrix object. If NULL, then fitted values are returned
#' @param ... not used
#' @return An object depending on the type argument
#' @importFrom bigmemory is.big.matrix as.big.matrix
#' @importFrom stats coef fitted pt printCoefmat
#' @export
#' @examples
#'
#' library(bigmemory)
#'
#' nrows <- 50000
#' ncols <- 50
#' bkFile <- "bigmat3.bk"
#' descFile <- "bigmatk3.desc"
#' bigmat <- filebacked.big.matrix(nrow=nrows, ncol=ncols, type="double",
#'                                 backingfile=bkFile, backingpath=".",
#'                                 descriptorfile=descFile,
#'                                 dimnames=c(NULL,NULL))
#'
#' # Each column value with be the column number multiplied by
#' # samples from a standard normal distribution.
#' set.seed(123)
#' for (i in 1:ncols) bigmat[,i] = rnorm(nrows)*i
#'
#' y <- rnorm(nrows) + bigmat[,1]
#'
#' system.time(lmr1 <- bigLm(bigmat, y))
#'
#'
#' preds <- predict(lmr1, newdata = bigmat)
#'
predict.bigLm <- function(object, newdata = NULL, ...) {
    if (is.null(newdata)) {
        y <- fitted(object)
    } else {
        stopifnot(is.big.matrix(newdata))
        y <- as.vector(newdata %*% coef(object))
    }
    y
}
