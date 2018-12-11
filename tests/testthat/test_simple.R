library(rstpm2)

context("Missing values")
brcancer2 <- rstpm2::brcancer
beta <- c(-7.25671730761628, -0.361401190022274, 4.76100002646889, 
          11.5713803414823, 4.56781356692313)
##
brcancer2$rectime[1] <- NA
fit1 <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer2)
## expect_equal(className(fit1), "stpm2")
expect_lt(max(abs(coef(fit1)-beta)), 1e-13)
expect_equal(length(predict(fit1)), 686)

brcancer2$rectime[1] <- -1
fit1 <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer2)
## expect_equal(className(fit1), "stpm2")
expect_lt(max(abs(coef(fit1)-beta)), 1e-13)
expect_equal(length(predict(fit1)), 686)

