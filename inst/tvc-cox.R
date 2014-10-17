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
lung1 <- lung3 <- within(subset(lung,!is.na(ph.ecog)), {
    event <- ifelse(status==1,1,0)
    x <- ifelse(ph.ecog==0,0,1)
})
lung2 <- lung1[rep(1:nrow(lung1),each=1000),]
lung3 <- lung3[rep(1:nrow(lung3),each=10),]
lung3$time <- lung3$time+rnorm(length(lung3$time),0,0.001)
lung3 <- lung3[with(lung3,order(time,-status)),]
## write.dta(lung1,file="~/work/lung1.dta")
## write.dta(lung2,file="~/work/lung2.dta")
## write.dta(lung3,file="~/work/lung3.dta")
system.time(fit0 <- coxph(Surv(time, event) ~ x , data=lung3))
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
system.time(fit <- coxph(Surv(time, event) ~ x + tt(x), data=lung3,
             tt=function(x,t,...) x*log(t)))
system.time(fit <- coxph(Surv(time, event) ~ x + tt(x), data=lung3,
             tt=function(x,t,...) x*log(t)))
(coxph(Surv(time, event) ~ x, data=lung3))

require(rstpm2)
fit0 <- pstpm2(Surv(time, event) ~ x, data=lung3, smooth.formula=~s(log(time)))
fit <- pstpm2(Surv(time, event) ~ 1, data=lung3, smooth.formula=~s(log(time))+s(log(time),by=x),sp.init=c(1,1))
fit2 <- pstpm2(Surv(time, event) ~ x, data=lung3, smooth.formula=~s(log(time))+x:log(time),sp.init=1)
summary(fit2)

i <- -(1:10)
coef1 <- coef(fit)[i]
statistic <- as.numeric(coef1 %*% solve(vcov(fit)[i,i]) %*% coef1)
1-pchisq(statistic,length(i))

fit0 <- stpm2(Surv(time, event) ~ x, data=lung3, df=3)
fit <- stpm2(Surv(time, event) ~ 1, data=lung3, tvc.formula=~x:ns(log(time),df=3))
fit2 <- stpm2(Surv(time, event) ~ x, data=lung3, tvc.formula=~x:log(time), df=4)
summary(fit2)
anova(fit2,fit0)

## Wald test by hand
i <- 5:7
coef1 <- coef(fit)[i]
statistic <- as.numeric(coef1 %*% solve(vcov(fit)[i,i]) %*% coef1)
1-pchisq(statistic,length(i))



system.time(fit <- pstpm2(Surv(time, event) ~ x, data=lung3, smooth.formula=~s(log(time))+log(time):x,sp.init=1))
library('Rcpp')
library('inline')
rcpp_inc <- '
using namespace Rcpp;
using namespace arma;
'
src <- '
    List largs = as<List>(args);
    vec time   = as<vec>(largs["time"]); // length n
    vec event  = as<vec>(largs["event"]); // length n
    vec x      = as<vec>(largs["x"]); // one covariate, length n
    vec beta   = as<vec>(largs["beta"]); // length 2
    int n      = time.size();
    double llike = 0.0, lsum;
    vec eta;
    for (int i=0; i<n; ++i) {
      if (event(i) == 1 && (i==0 || time(i) != time(i-1))) {
	eta = beta(0)*x(span(i,n-1)) + beta(1)*log(time(i))*x(span(i,n-1));
        lsum = log(sum(exp(eta)));
        for (int j=i; j<n && time(j)==time(i) && event(j)==1; ++j) {
          llike += eta(j-i) - lsum;
        }
      }
    }
    return wrap(llike);
'
fn <- cxxfunction(signature(args="list"), src, plugin='RcppArmadillo', rcpp_inc)
src <- '
    List largs = as<List>(args);
    vec time   = as<vec>(largs["time"]); // length n
    vec event  = as<vec>(largs["event"]); // length n
    vec x      = as<vec>(largs["x"]); // one covariate, length n
    vec beta   = as<vec>(largs["beta"]); // length 2
    int n      = time.size();
    vec grad(beta.size());
    vec risk;
    vec lsum(2);
    for (int i=0; i<n; ++i) {
      if (event(i) == 1 && (i==0 || time(i) != time(i-1))) {
	risk = exp(beta(0)*x(span(i,n-1)) + beta(1)*log(time(i))*x(span(i,n-1)));
        lsum(0) = sum(x(span(i,n-1)) % risk)/sum(risk);
        lsum(1) = sum(log(time(i))*x(span(i,n-1)) % risk)/sum(risk);
        for (int j=i; j<n && time(j)==time(i) && event(j)==1; ++j) {
	  grad(0) += x(j) - lsum(0);
	  grad(1) += x(j)*log(time(i)) - lsum(1);
        }
      }
    }
    return wrap(grad);
'
gr <- cxxfunction(signature(args="list"), src, plugin='RcppArmadillo', rcpp_inc)


src <- '
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
double fn(NumericVector time_, NumericVector event_, NumericVector x_, NumericVector beta_) { 
    vec event = as<vec>(wrap(event_));
    vec time = as<vec>(wrap(time_));
    vec x = as<vec>(wrap(x_));
    vec beta = as<vec>(wrap(beta_));
    int n      = time.size();
    double llike = 0.0, lsum;
    vec eta;
    for (int i=0; i<n; ++i) {
      if (event(i) == 1 && (i==0 || time(i) != time(i-1))) {
	eta = beta(0)*x(span(i,n-1)) + beta(1)*log(time(i))*x(span(i,n-1));
        lsum = log(sum(exp(eta)));
        for (int j=i; j<n && time(j)==time(i) && event(j)==1; ++j) {
          llike += eta(j-i) - lsum;
        }
      }
    }
    return llike;
}'
sourceCpp(code=src)
with(lung1,fn(time,event,x,beta))


lung1 <- lung1[with(lung1,order(time,-event)),]
fit0 <- coxph(Surv(time, event) ~ x , data=lung1)
beta <- c(coef(fit0),0)
(fit <- coxph(Surv(time, event) ~ x + tt(x), data=lung1,
              tt=function(x,t,...) x*log(t)))
args <- with(lung1,list(time=time,event=event,x=x,beta=beta))
fn(args)
gr(args)
.Call("test_cox_tvc2",args,package="rstpm2")

system.time(fit3 <- optim(par=beta,
                          fn = function(beta) -fn(list(time=lung1$time,event=lung1$event,x=lung1$x,beta=beta)),
                          hessian=TRUE, control=list(abstol=1e-8,reltol=1e-8,trace=1)))
fit3$par
system.time(fit4 <- optim(par=beta,
                          fn = function(beta) -fn(list(time=lung1$time,event=lung1$event,x=lung1$x,beta=beta)),
                          gr = function(beta) -gr(list(time=lung1$time,event=lung1$event,x=lung1$x,beta=beta)),
                          method="BFGS",
                          hessian=TRUE, control=list(abstol=1e-8,reltol=1e-8,trace=1)))
fit4$par

system.time(fit3 <- optim(par=beta,
                          fn = function(beta) -fn(list(time=lung3$time,event=lung3$event,x=lung3$x,beta=beta)),
                          hessian=TRUE, control=list(abstol=1e-8,reltol=1e-8,trace=1)))
fit3$par

lung2 <- lung2[with(lung2,order(time,-event)),]
args <- with(lung2,list(time=time,event=event,x=x,beta=beta))
fn(args)
system.time(fit3 <- optim(par=beta,
                          fn = function(beta) -fn(list(time=lung2$time,event=lung2$event,x=lung2$x,beta=beta)),
                          hessian=TRUE, control=list(abstol=1e-8,reltol=1e-8,trace=1)))
fit3$par



library('Rcpp')
library('inline')
rcpp_inc <- '
using namespace Rcpp;
using namespace arma;
'
src <- '
ivec index = linspace<ivec>(0,10,11);
return(wrap(find(index>5 && index/8==1)));
'
(cxxfunction(signature(), src, plugin='RcppArmadillo', rcpp_inc))()
src <- '
    List largs = as<List>(args);
    vec time   = as<vec>(largs["time"]); // length n
    vec event  = as<vec>(largs["event"]); // length n
    vec x      = as<vec>(largs["x"]); // one covariate, length n
    vec beta   = as<vec>(largs["beta"]); // length 2
    int n      = time.size();
    double llike = 0.0;
    vec eta;
    for (int i=0; i<n; ++i) {
      if (event(i) == 1) {
	eta = beta(0)*x(span(i,n-1)) + beta(1)*log(time(i))*x(span(i,n-1));
	llike += eta(0) - log(sum(exp(eta)));
      }
    }
    return wrap(llike);
'
fn <- cxxfunction(signature(args="list"), src, plugin='RcppArmadillo', rcpp_inc)


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
