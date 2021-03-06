% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fast_big_lm.R
\name{summary.bigLm}
\alias{summary.bigLm}
\title{summary method for bigLm fitted objects}
\usage{
\method{summary}{bigLm}(object, ...)
}
\arguments{
\item{object}{bigLm fitted object}

\item{...}{not used}
}
\value{
a summary.bigLm object
}
\description{
summary method for bigLm fitted objects
}
\examples{

library(bigmemory)

nrows <- 50000
ncols <- 15
bkFile <- "bigmat4.bk"
descFile <- "bigmatk4.desc"
bigmat <- filebacked.big.matrix(nrow=nrows, ncol=ncols, type="double",
                                backingfile=bkFile, backingpath=".",
                                descriptorfile=descFile,
                                dimnames=c(NULL,NULL))

# Each column value with be the column number multiplied by
# samples from a standard normal distribution.
set.seed(123)
for (i in 1:ncols) bigmat[,i] = rnorm(nrows)*i

y <- rnorm(nrows) + bigmat[,1] - bigmat[,2]

system.time(lmr1 <- bigLm(bigmat, y))

summary(lmr1)


}

