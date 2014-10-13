## Times are shown for my computer
require(survival)
require(foreign)
lung1 <- subset(survival::lung, !is.na(ph.ecog)) # one case with missing ph.ecog
lung2 <- lung1[rep(1:nrow(lung1),each=10),] # n=227000 - and many ties
## write.dta(lung1,file="~/work/lung1.dta")
## write.dta(lung2,file="~/work/lung2.dta")
system.time(fit0 <- coxph(Surv(time, status) ~ ph.ecog, data=lung2)) # 0.5s
system.time(fit <- coxph(Surv(time, status) ~ ph.ecog + tt(ph.ecog), data=lung1,
             tt=function(x,t,...) x*log(t))) # 300s - SLOW
anova(fit)
system.time(cox.zph(fit0)) # 2.2s

refresh
require(survival)
require(rstpm2)
require(foreign)
set.seed(12345)
lung3 <- within(subset(lung,!is.na(ph.ecog)), {
    event <- ifelse(status==1,1,0)
    x <- ifelse(ph.ecog==0,0,1)
})
lung3 <- lung3[rep(1:nrow(lung3),each=10),]
lung3$time <- lung3$time+rnorm(length(lung3$time),0,0.001)
lung3 <- lung3[with(lung3,order(time,-status)),]
## write.dta(lung3,file="~/work/lung3.dta")
system.time(fit0 <- coxph(Surv(time, event) ~ x , data=lung3))
## system.time(fit <- coxph(Surv(time, event) ~ x + tt(x), data=lung3,
##                          tt=function(x,t,...) x*log(t)))
beta <- c(coef(fit0),0)
## .Call("test_cox_tvc",with(lung3,list(time=time,event=event,x=x,beta=beta)),package="rstpm2")
## .Call("test_cox_tvc2",with(lung3,list(time=time,event=event,x=x,beta=beta)),package="rstpm2")
## .Call("test_cox_tvc2_grad",with(lung3,list(time=time,event=event,x=x,beta=beta)),package="rstpm2")

system.time(fit3 <- optim(par=beta,
                          fn = function(beta) -.Call("test_cox_tvc2",
                              list(time=lung3$time,event=lung3$event,x=lung3$x,beta=beta),package="rstpm2"),
                          gr = function(beta) -.Call("test_cox_tvc2_grad",
                              list(time=lung3$time,event=lung3$event,x=lung3$x,beta=beta),package="rstpm2"),
                          method="BFGS", hessian=TRUE, control=list(trace=1)))
fit3$par
##
system.time(fit4 <-
            .Call("test_cox_tvc3",
                  list(time=lung3$time,event=lung3$event,x=lung3$x,beta=beta),package="rstpm2"))
fit4


## require(devtools)
## install_github("mclements/rstpm2", ref="develop", quick=TRUE)
require(rstpm2)
system.time(fit0 <- stpm2(Surv(time, status) ~ ph.ecog, data=lung2)) # 6.5s
system.time(fit <- stpm2(Surv(time, status) ~ ph.ecog, data=lung2, tvc=list(ph.ecog=3))) # 9.3s
anova(fit,fit0)
plot(fit,newdata=data.frame(ph.ecog=0),type="hazard")
plot(fit,newdata=data.frame(ph.ecog=1),type="hazard",add=TRUE,line.col="blue",ci=TRUE)
legend("topleft", legend=c("ph.ecog=0","ph.ecog=1"),lty=1,col=c("black","blue"),bty="n")

system.time(fit <- pstpm2(Surv(time, status) ~ 1, data=lung2, smooth.formula=~s(log(time))+s(log(time),by=ph.ecog),sp.init=c(1,100))) # 400s
plot(fit,newdata=data.frame(ph.ecog=0),type="hazard",ylim=c(0,0.005))
plot(fit,newdata=data.frame(ph.ecog=1),type="hazard",add=TRUE,line.col="blue",ci=TRUE)

## install.packages("CoxRidge")
require(CoxRidge)
Ft <- with(lung2, cbind(1,bs(time)))
X <- subset(lung2, select=ph.ecog)
system.time(fit.dr <- Dynamic.Ridge(lung2$time,lung2$status,X,Ft=Ft,lambda=100,fun="simple")) # FAILS
fit.dr
