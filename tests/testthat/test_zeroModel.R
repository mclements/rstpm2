library(rstpm2)

expect_eps <- function(expr, value, eps=1e-7)
    expect_lt(max(abs(expr-value)),eps)

context("zeroModel")
##
test_that("base", {
    x <- 1:10
    y <- c(1:9,11)
    d <- data.frame(x,y)
    fit <- zeroModel(lm(y~x,data=d))
    expect_eps(coef(fit), c(0,0), 1e-10)
    expect_eps(vcov(fit), matrix(0,2,2), 1e-10)
    ## expect_eps(predict(fit,newdata=d), rep(0,10), 1e-10) # zeroModel class not exported
})

context("hrModel")
##
test_that("base", {
    x <- 1:10
    y <- c(1:9,11)
    fit <- hrModel(glm(y~x,family=poisson),2,ci=c(1,4))
    expect_eps(coef(fit), c(0.4577646, 0.2007416, 0.6931472), 1e-5)
    expect_eps(vcov(fit),
               matrix(c(0.148392636125595, -0.0185062979900643, 0,
                        -0.0185062979900643, 0.00262367771332864, 0,
                        0, 0, 0.125070457954665),3,3), 1e-10)
    expect_eps(predict(fit), predict.glm(fit$base,type="haz")*2, 1e-10)
    expect_eps(predict(fit,type="gradh"),
               cbind(predict.glm(fit$base,type="gradh")*2,
                     predict.glm(fit$base,type="haz")*2),
               1e-10)
})
