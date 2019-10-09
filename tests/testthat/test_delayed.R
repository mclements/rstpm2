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

