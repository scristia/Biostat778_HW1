##########################################
## A faster linear regression.
## Arguments:
##      X = an nxp matrix
##      y = a vector of length n
##      na.rm
##########################################
fastlm <- function(X, y, na.rm=FALSE) {
    ## dumb way
    ## take advantage of X'X being symmetric
    n <- nrow(X)
    p <- ncol(X)
    U <- chol(t(X)%*%X)

    b <- backsolve(U, forwardsolve(t(U), t(X)%*%y))
    return(b)
    #sigma <- 1/(n-p) * (t(y) %*% y - t(b) %*% t(X) %*% y)
    #return(list("coefficients" = b, "vcov" = sigma))
}
fastlm.svd <- function(X, y, na.rm=FALSE) {
    ## dumb way
    ## take advantage of X'X being symmetric
    m <- svd(t(X)%*%X)
    b <- m$v %*% diag(m$d^(-1)) %*% t(m$u) %*% t(X) %*% y

    return(b)
}
fastlm.bad <- function(X, y, na.rm=FALSE) {
    ## dumb way
    ## take advantage of X'X being symmetric
    b <- solve(t(X)%*%X, t(X)%*%y)

    return(b)
}
system.time(fastlm(X,y))
system.time(fastlm.svd(X,y))
system.time(fastlm.bad(X,y))
str(fastlm(X,y))
ans <- solve(t(X)%*%X, t(X)%*%y)

set.seed(2)
## Generate predictor matrix
n <- 100000
p <- 500
X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))

## Coefficents
b <- rnorm(p)

## Response
y <- X %*% b + rnorm(n)
