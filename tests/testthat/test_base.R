library(rstpm2)

## for coping with weird test behaviour from CRAN and R-devel
.CRAN <- FALSE
slow <- FALSE

expect_eps <- function(expr, value, eps=1e-7)
    expect_lt(max(abs(expr-value)),eps)

context("stpm2")
##
test_that("base", {
    fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer)
    expect_eps(coef(fit)[2], -0.361403, 1e-5)
    fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,
                 control=list(optimiser="NelderMead"))
    expect_eps(coef(fit)[2], -0.3608321, 1e-5)
})
test_that("Comparison with Stata", {
    fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,
                 smooth.formula=~nsx(log(rectime),df=3,stata=TRUE))
    expect_eps(coef(fit)[2], -0.3614357, 1e-5)
})
test_that("Cure", {
    fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,
                 smooth.formula=~nsx(log(rectime),df=3,cure=TRUE))
    expect_eps(coef(fit)[2], -0.3564268, 1e-5)
    fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,cure=TRUE)
    expect_eps(coef(fit)[2], -0.3564268, 1e-5)
})

context("stpm2 + frailty")
##
test_that("base", {
    set.seed(12345)
    id <- rep(1:200,each=2)
    x <- rnorm(400)
    mu.cluster <- rnorm(200)
    y <- rexp(400, exp(mu.cluster[id]+x))
    d <- data.frame(y,e=1,x,id)
    fit <- stpm2(Surv(y,e)~x+cluster(id),data=d)
    expect_eps(coef(fit)[2], 0.9256663, 1e-5)
    expect_eps(coef(fit)[6], -0.4698344, 1e-5)
    fit <- stpm2(Surv(y,e)~x+cluster(id),data=d,RandDist="LogN")
    expect_eps(coef(fit)[2], 0.93819784, 1e-5)
    expect_eps(coef(fit)[6], 0.02906888, 1e-5)
})

context("pstpm2")
##
test_that("base", {
    fit <- pstpm2(Surv(rectime,censrec==1)~hormon,data=brcancer)
    expect_eps(coef(fit)[2], -0.3650564, 1e-5)
    fit <- pstpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,
                  control=list(optimiser="NelderMead"))
    expect_eps(coef(fit)[2], -0.3652168, 1e-5)
})

if (slow) {
    context("pstpm2 + frailty")
    ##
    test_that("base", {
        set.seed(12345)
        id <- rep(1:200,each=2)
        x <- rnorm(400)
        mu.cluster <- rnorm(200)
        y <- rexp(400, exp(mu.cluster[id]+x))
        d <- data.frame(y,e=1,x,id)
        fit <- pstpm2(Surv(y,e)~x+cluster(id),data=d)
        expect_eps(coef(fit)[2], 0.9071705, 1e-5)
        expect_eps(coef(fit)[12], -5.303905e-01, 1e-5)
        fit <- pstpm2(Surv(y,e)~x+cluster(id),data=d,RandDist="LogN")
        expect_eps(coef(fit)[2], 9.334020e-01, 1e-5)
        expect_eps(coef(fit)[6], 1.426793e-02, 1e-5)
    })
}
