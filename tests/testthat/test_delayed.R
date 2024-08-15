library(rstpm2)

## for coping with weird test behaviour from CRAN and R-devel
.CRAN <- FALSE

expect_eps <- function(expr, value, eps=1e-7)
    expect_lt(max(abs(expr-value)),eps)

context("Delayed entry - stpm2")
##
test_that("Comparison with Stata", {
    brcancer2 <- transform(brcancer,startTime=ifelse(hormon==0,rectime/2,0))
    fit <- stpm2(Surv(startTime,rectime,censrec==1)~hormon,data=brcancer2,
                 smooth.formula=~nsx(log(rectime),df=3,stata=TRUE))
    expect_eps(coef(fit)[2], -1.162504, 1e-5)
})
test_that("pstpm2", {
    brcancer2 <- transform(brcancer,startTime=ifelse(hormon==0,rectime/2,0))
    fit <- pstpm2(Surv(startTime,rectime,censrec==1)~hormon,data=brcancer2)
    expect_eps(coef(fit)[2], -1.1881852, 1e-5)
    ## fit <- pstpm2(Surv(startTime,rectime,censrec==1)~hormon,data=brcancer2,
    ##               control=list(optimiser="NelderMead"))
    ## expect_eps(coef(fit)[2], -1.193484, 1e-5)
})

context("Delayed entry - aft")
##
test_that("All values zero or one", {
    brcancer2 <- transform(rstpm2::brcancer,startTime=0)
    fit0 <- aft(Surv(rectime,censrec==1)~hormon,data=brcancer2)
    fit1 <- aft(Surv(startTime,rectime,censrec==1)~hormon,data=brcancer2)
    expect_true(all(coef(fit0) == coef(fit1)))
    brcancer2 <- transform(rstpm2::brcancer,startTime=1)
    fit2 <- aft(Surv(startTime,rectime,censrec==1)~hormon,data=brcancer2)
    diff <- abs(coef(fit2)-coef(fit1))
    expect_true(min(diff)>1e-7)
    expect_true(max(diff)<1e-5)
    set.seed(12345)
    brcancer3 <- transform(rstpm2::brcancer,startTime=rbinom(nrow(rstpm2::brcancer),1,0.5))
    fit3 <- aft(Surv(startTime,rectime,censrec==1)~hormon,data=brcancer3)
    diff <- abs(coef(fit2)-coef(fit1))
    expect_true(min(diff)>1e-7)
    expect_true(max(diff)<1e-5)
})
