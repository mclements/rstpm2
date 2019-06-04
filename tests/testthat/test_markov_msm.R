library(rstpm2)

## for coping with weird test behaviour from CRAN and R-devel
.CRAN <- FALSE
## pstpm2+frailty models are slow
slow <- FALSE

expect_eps <- function(expr, value, eps=1e-7)
    expect_lt(max(abs(expr-value)),eps)

context("markov_msm")
##
if (slow) {
    test_that("gam", {
        library(mgcv)
        library(survival)
        set.seed(12345)
        t <- rweibull(1e3,shape=2)
        d <- data.frame(t,e=1)
        tsplit <- survival::survSplit(Surv(t,e)~1, data=d, cut=seq(0,6,length=100))
        tsplit <- transform(tsplit, dt=t-tstart)
        fit <- gam(e~s(log(t))+offset(log(dt)), data=tsplit, family=poisson)
        fit <- markov_msm(list(fit), trans=matrix(c(NA,1,NA,NA),2,2,TRUE), newdata=data.frame(dt=1),
                          t=c(0,2),tmvar="t")
        expect_eps(as.data.frame(fit)$P[4], 0.9786583, 1e-5)
    })
}
test_that("glm", {
    library(survival)
    set.seed(12345)
    t <- rweibull(1e3,shape=2)
    d <- data.frame(t,e=1)
    tsplit <- survival::survSplit(Surv(t,e)~1, data=d, cut=seq(0,6,length=100))
    tsplit <- transform(tsplit, dt=t-tstart)
    fit <- glm(e~ns(log(t),df=3)+offset(log(dt)), data=tsplit, family=poisson)
    fit <- markov_msm(list(fit), trans=matrix(c(NA,1,NA,NA),2,2,TRUE), newdata=data.frame(dt=1),
                      t=c(0,2),tmvar="t")
    expect_eps(as.data.frame(fit)$P[4], 0.977278, 1e-5)
})

