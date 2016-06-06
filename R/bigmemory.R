

big.crossprod <- function(x)
{
  if (!is.big.matrix(x))
  {
    stop("object must be a big.matrix object")
  }
  .Call("crossprod_big", x@address)
}

big.colSums <- function(x)
{
  if (!is.big.matrix(x))
  {
    stop("object must be a big.matrix object")
  }
  .Call("colsums_big", x@address)
}

big.colMax <- function(x)
{
  if (!is.big.matrix(x))
  {
    stop("object must be a big.matrix object")
  }
  .Call("colmax_big", x@address)
}

big.colMin <- function(x)
{
  if (!is.big.matrix(x))
  {
    stop("object must be a big.matrix object")
  }
  .Call("colmin_big", x@address)
}


setMethod("%*%",signature(x="big.matrix", y="vector"),
          function(x, y) 
          {
            if(dim(x)[2] != length(y)) stop("non-conformant matrices")
            return( .Call("prod_vec_big", x@address, y) )
          },
          valueClass="vector"
)

setMethod("%*%",signature(x="vector", y="big.matrix"),
          function(x, y) 
          {
            if(dim(y)[2] != length(x)) stop("non-conformant matrices")
            return( .Call("prod_vec_big_right", x, y@address) )
          },
          valueClass="vector"
)


