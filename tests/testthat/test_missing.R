library(rstpm2)

## context("Missing values")
expect_eps <- function(expr, value, eps=1e-7)
    expect_lt(max(abs(expr-value)),eps)

context("Missing event time")
brcancer2 <- rstpm2::brcancer
beta1 <- c(-7.24333895936283, -0.359488397062925, 4.75357528198452, 
           11.5332943620568, 4.5659764943589)
##
brcancer2$rectime[1] <- NA
expect_warning(fit1 <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer2),
               "Some event times are NA")
expect_eps(coef(fit1),beta1, 1e-13)
expect_length(predict(fit1), 685)

context("Invalid event time")
brcancer2$rectime[1] <- -1
expect_warning(fit1 <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer2),
               "Some event times <= 0")
expect_eps(coef(fit1),beta1, 1e-13)
expect_length(predict(fit1), 685)

context("Missing covariate")
brcancer2 <- rstpm2::brcancer
brcancer2$hormon[1] <- NA
fit1 <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer2)
expect_eps(coef(fit1),beta1, 1e-13)
expect_length(predict(fit1), 685)

context("Missing weight")
brcancer2 <- transform(rstpm2::brcancer, w=1)
brcancer2$w[1] <- NA
fit1 <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer2,weights=w)
expect_eps(coef(fit1),beta1, 1e-13)
expect_length(predict(fit1), 685)

context("Predictions with missing values")
test <- c(surv=0.7230207, fail=1-0.7230207, haz=0.000335262)
for(name in names(test))
    expect_eps(predict(fit1,newdata=data.frame(hormon=1,rectime=1000),type=name),
               test[name], 1e-6)
test <- c(hr=0.6980334, hdiff=-0.000101238)
for(name in names(test))
    expect_eps(predict(fit1,newdata=data.frame(hormon=1,rectime=1000),type=name,var="hormon"),
               test[name], 1e-6)
