library(rstpm2)

## context("Missing values")

context("Missing event time")
brcancer2 <- rstpm2::brcancer
beta1 <- c(-7.24333895936283, -0.359488397062925, 4.75357528198452, 
           11.5332943620568, 4.5659764943589)
##
brcancer2$rectime[1] <- NA
expect_warning(fit1 <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer2),
               "Some event times are NA")
expect_lt(max(abs(coef(fit1)-beta1)), 1e-13)
expect_equal(length(predict(fit1)), 685)

context("Invalid event time")
brcancer2$rectime[1] <- -1
expect_warning(fit1 <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer2),
               "Some event times <= 0")
expect_lt(max(abs(coef(fit1)-beta1)), 1e-13)
expect_equal(length(predict(fit1)), 685)

context("Missing covariate")
brcancer2 <- rstpm2::brcancer
brcancer2$hormon[1] <- NA
fit1 <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer2)
expect_lt(max(abs(coef(fit1)-beta1)), 1e-13)
expect_equal(length(predict(fit1)), 685)

context("Missing weight")
brcancer2 <- transform(rstpm2::brcancer, w=1)
brcancer2$w[1] <- NA
fit1 <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer2,weights=w)
expect_lt(max(abs(coef(fit1)-beta1)), 1e-13)
expect_equal(length(predict(fit1)), 685)

context("Predictions with missing values")
test <- c(surv=0.7230207, fail=1-0.7230207, haz=0.000335262)
for(this in names(test))
    expect_lt(abs(predict(fit1,newdata=data.frame(hormon=1,rectime=1000),type=this)-test[this]),
              1e-6)
