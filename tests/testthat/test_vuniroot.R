library(rstpm2)

## for coping with weird test behaviour from CRAN and R-devel
.CRAN <- FALSE
## pstpm2+frailty models are slow
slow <- FALSE

expect_eps <- function(expr, value, eps=1e-7)
    expect_lt(max(abs(expr-value)),eps)

context("vuniroot")
##
test_that("examples", {
    ## some platforms hit zero exactly on the first step:
    ## if so the estimated precision is 2/3.
    f <- function (x, a) x - a
    xmin <- vuniroot(f, lower=c(0, 0), upper=c(1,1), tol = 0.0001, a = c(1/3,2/3))
    expect_eps(xmin$root, c(1/3,2/3), 1e-10)
    ## handheld calculator example: fixed point of cos(.):
    expect_eps(vuniroot(function(x) cos(x) - x, lower = -pi, upper = pi, tol = 1e-9)$f.root,
               0, 1e-10)
    expect_eps(vuniroot(function(x) x*(x^2-1) + .5, lower = -2, upper = 2,
            tol = 0.0001)$f.root, 0, 1e-5)
    expect_eps(vuniroot(function(x) x*(x^2-1) + .5, lower = -2, upper = 2,
                        tol = 1e-10)$f.root, 0, 1e-10)
    ## Find the smallest value x for which exp(x) > 0 (numerically):
    expect_eps((r <- vuniroot(function(x) 1e80*exp(x) - 1e-300, cbind(-1000, 0), tol = 1e-15))$root,
               -745.133219101941, 1e-5)
    ##--- vuniroot() with new interval extension + checking features: --------------
    f1 <- function(x) (121 - x^2)/(x^2+1)
    f2 <- function(x) exp(-x)*(x - 12)
    expect_error(vuniroot(f1, cbind(0,10)),
                 "f[(][)] values at end points not of opposite sign")
    expect_error(vuniroot(f1, cbind(0,2)),
                 "f[(][)] values at end points not of opposite sign")
    ## where as  'extendInt="yes"'  simply first enlarges the search interval:
    u1 <- vuniroot(f1, cbind(0,10),extendInt="yes")
    u2 <- vuniroot(f2, cbind(0,2), extendInt="yes")
    expect_eps(u1$root, 11, 1e-4)
    expect_eps(u2$root, 12, 1e-5)
    ## The *danger* of interval extension:
    ## No way to find a zero of a positive function, but
    ## numerically, f(-|M|) becomes zero :
    expect_error(vuniroot(exp, cbind(0,2), extendInt="yes"),
                 "did not succeed extending the interval endpoints for f[(]lower[)] [*] f[(]upper[)] <= 0")
    ## Nonsense example (must give an error):
    expect_error(vuniroot(function(x) 1, cbind(0,1), extendInt="yes"),
                 "no sign change found in 1000 iterations")
    ## Convergence checking :
    sinc <- function(x) ifelse(x == 0, 1, sin(x)/x)
    expect_warning(vuniroot(sinc, cbind(0,5), extendInt="yes", maxiter=4),
                   "_NOT_ converged in 4 iterations")
    ## now with  check.conv=TRUE, must signal a convergence error :
    expect_error(vuniroot(sinc, cbind(0,5), extendInt="yes", maxiter=4, check.conv=TRUE),
                   "_NOT_ converged in 4 iterations")
    ## Weibull cumulative hazard (example origin, Ravi Varadhan):
    cumhaz <- function(t, a, b) b * (t/b)^a
    froot <- function(x, u, a, b) cumhaz(x, a, b) - u
    set.seed(12345)
    n <- 10
    u <- -log(runif(n))
    a <- 1/2
    b <- 1
    ## Find failure times
    ru <- vuniroot(froot, u=u, a=a, b=b, interval= cbind(rep(1.e-14,n), rep(1e4,n)),
                   extendInt="yes")$root
    ru2 <- vuniroot(froot, u=u, a=a, b=b, interval= cbind(rep(0.01,n), rep(10,n)),
                    extendInt="yes")$root
    expect_eps(ru, ru2, 1e-4)
    r1 <- vuniroot(froot, u= 0.99, a=a, b=b, interval= cbind(0.01, 10),
                   extendInt="up")
    expect_eps(0.99, cumhaz(r1$root, a=a, b=b), 1e-8)
    expect_error(vuniroot(froot, u= 0.99, a=a, b=b, interval= cbind(0.1, 10), extendInt="down"),
                 "no sign change found in 1000 iterations")
})
