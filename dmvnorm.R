dmvnorm <- function(x, mu, S, log = TRUE) {
    ## x = data vector, mu = mean, S = covariance matrix
    ## if x numeric array, make row vector
    if(!is.matrix(x)) x <- matrix(x, nrow=1, ncol=ncol(S))
    k <- ncol(S)
    ## need: fast quadratic form
    ## fast determinant of positive definite matrix
    U <- try(chol(S), silent=TRUE)
    if(class(U) == "try-error") stop("S is not positive definite")

    d <- diag(U)
    ## check if positive definite
    logd <- sum(log(d))
    b <- crossprod(forwardsolve(U, t(x-mu), upper.tri=TRUE))
#    b <- crossprod(backsolve(U, t(x-mu)))
    logf  <- -k/2 * log(2*pi) - logd - 0.5 * b
    if(log) return(logf)
    else return(exp(logf))
}

is.pos.def <- function(x) {
    ## check symmetric
    return(sum(x == t(x)) == (nrow(x)^2))
}
