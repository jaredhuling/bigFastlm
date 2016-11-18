
## bigFastlm [![version](http://www.r-pkg.org/badges/version/bigFastlm)](https://cran.r-project.org/package=bigFastlm)
A reimplementation of the fastLm functionality of RcppEigen for
big.matrix objects for fast out-of-memory linear model fitting


### Build Status
|  OS                   | Build           |
|-----------------------|-----------------|
| Linux x86_64          | [![Build Status](https://travis-ci.org/jaredhuling/bigFastlm.svg?branch=master)](https://travis-ci.org/jaredhuling/bigFastlm)      | 
| Windows x86_64        | [![Appveyor Build Status](https://ci.appveyor.com/api/projects/status/github/jaredhuling/bigFastlm?branch=master&svg=true)](https://ci.appveyor.com/project/jaredhuling/bigFastlm)     |



### Installation

Install using the **devtools** package (RcppEigen and bigmemory must be installed first as well):

```r
devtools::install_github("jaredhuling/bigFastlm")
```


### Usage

```r
library(bigFastlm)
library(bigmemory)

nrows <- 50000
ncols <- 50
bkFile <- "bigmat.bk"
descFile <- "bigmatk.desc"
bigmat <- filebacked.big.matrix(nrow=nrows, ncol=ncols, type="double",
                                backingfile=bkFile, backingpath=".",
                                descriptorfile=descFile,
                                dimnames=c(NULL,NULL))

set.seed(123)
for (i in 1:ncols) bigmat[,i] = rnorm(nrows)*i

y <- rnorm(nrows) + bigmat[,1]

system.time(lmr1 <- bigLm(bigmat, y))

summary(lmr1)

predictions <- predict(lmr1, newdata = bigmat)
```
