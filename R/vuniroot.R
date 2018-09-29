vuniroot <- 
function (f, ..., lower, upper, 
    f.lower = f(lower, ...), f.upper = f(upper, ...), check.conv = FALSE, tol = .Machine$double.eps^0.25, 
    maxiter = 1000, trace = 0) 
{
    stopifnot(all(is.numeric(lower), is.numeric(upper), is.numeric(f.lower), is.numeric(f.upper)))
    stopifnot(all.equal(length(lower), length(upper), length(f.lower), length(f.upper)))
    stopifnot(!any(is.na(lower)))
    stopifnot(!any(is.na(upper)))
    stopifnot(!any(is.na(f.lower)))
    stopifnot(!any(is.na(f.upper)))
    ## re-order lower and upper if lower>upper - too automagical?
    if (any(index <- (lower>upper))) {
        temp <- lower[index]
        lower[index] <- upper[index]
        upper[index] <- temp
        temp <- f.lower[index]
        f.lower[index] <- f.upper[index]
        f.upper[index] <- temp
    }
    maxiter <- as.integer(maxiter)
    fun <- function(x) f(x, ...)
    if (check.conv) {
        val <- tryCatch(.Call("vunirootRcpp", fun, lower, upper, maxiter, tol, PACKAGE="rstpm2"), 
            warning = function(w) w)
        if (inherits(val, "warning")) 
            stop("convergence problem in zero finding: ", conditionMessage(val))
    }
    else {
        val <- .Call("vunirootRcpp", fun, lower, upper, maxiter, tol, PACKAGE="rstpm2")
    }
    iter <- as.integer(val[[2L]])
    if (any(iter < 0)) {
        (if (check.conv) 
            stop
        else warning)(sprintf(ngettext(maxiter, "_NOT_ converged in %d iteration", 
            "_NOT_ converged in %d iterations"), maxiter), domain = NA)
        iter <- maxiter
    }
    list(root = val[[1L]], f.root = f(val[[1L]], ...), iter = iter) 
        #init.it = it), estim.prec = val[[3L]])
}
