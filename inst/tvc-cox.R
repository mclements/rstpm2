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
library('Rcpp')
library('inline')
rcpp_inc <- '
using namespace Rcpp;
using namespace arma;
'
set.seed(12345)
lung1 <- lung3 <- within(subset(lung,!is.na(ph.ecog)), {
    event <- ifelse(status==1,1,0)
    x <- ifelse(ph.ecog==0,0,1)
})
lung2 <- lung1[rep(1:nrow(lung1),each=1000),]
lung3 <- lung3[rep(1:nrow(lung3),each=10),]
lung3$time <- lung3$time+rnorm(length(lung3$time),0,0.001)
lung3 <- lung3[with(lung3,order(time,-status)),]
if (FALSE) {
    write.dta(lung1,file="~/work/lung1.dta")
    write.dta(lung2,file="~/work/lung2.dta")
    write.dta(lung3,file="~/work/lung3.dta")
}
system.time(fit0 <- coxph(Surv(time, event) ~ x , data=lung3))
beta <- betaInit <- c(coef(fit0),0)

tvc.coxph <- function(obj,var,method="logt") {
    stopifnot(attr(obj$y,"type") == "right")
    stopifnot(method == "logt")
    y <- as.matrix(obj$y)
    time <- y[,1]
    status <- y[,2]
    X <- as.matrix(model.matrix(obj))
    index <- order(time,-status)
    X <- X[index, , drop = FALSE]
    time <- time[index]
    status <- status[index]
    k <- match(attr(obj$terms,"term.labels"),var)
    beta <- c(coef(obj),0)
    names(beta) <- c(names(coef(obj)),sprintf("%s:log(t)",var))
    minuslogl <- function(beta) -.Call("test_cox_tvc2",
                                   list(time=time,event=status,X=X,beta=beta,k=k-1),package="rstpm2")
    gr <- function(beta) -.Call("test_cox_tvc2_grad",
                                list(time=time,event=status,X=X,beta=beta,k=k-1),package="rstpm2")
    parnames(minuslogl) <- parnames(gr) <- names(beta)
    fit <- mle2(start=beta,
                 minuslogl = minuslogl,
                 gr = gr,
                 method="BFGS", hessian=TRUE)
    fit@data <- list(object=obj)
    ## plot the result
    betak <- beta*0
    index2 <- c(k,length(betak))
    betak[index2] <- coef(fit)[index2]
    Xk <- cbind(X,log(time))
    Xk[,-index2] <- 0
    Xk[,k] <- 1
    fitted <- as.vector(Xk %*% betak)
    gd <- t(Xk)
    se.fit <- sqrt(colSums(gd * (vcov(fit) %*% gd)))
    matplot(time,exp(fitted+cbind(0,-1.96*se.fit,1.96*se.fit)),type="l",lty=c(1,2,2),col=1,ylab="Effect",log="y")
    abline(h=exp(coef(obj)[k]),lty=3)
    fit
}

summary(tvc.coxph(fit0,"x"))
    
system.time(fit3 <- optim(par=beta,
                          fn = function(beta) -.Call("test_cox_tvc2",
                              list(time=lung3$time,event=lung3$event,X=X,beta=beta,k=0),package="rstpm2"),
                          gr = function(beta) -.Call("test_cox_tvc2_grad",
                              list(time=lung3$time,event=lung3$event,X=X,beta=beta,k=0),package="rstpm2"),
                          method="BFGS", hessian=TRUE, control=list(trace=1)))


system.time(fit3 <- optim(par=beta,
                          fn = function(beta) -.Call("test_cox_tvc2",
                              list(time=lung3$time,event=lung3$event,X=cbind(lung3$x),beta=beta,k=0),package="rstpm2"),
                          gr = function(beta) -.Call("test_cox_tvc2_grad",
                              list(time=lung3$time,event=lung3$event,X=cbind(lung3$x),beta=beta,k=0),package="rstpm2"),
                          method="BFGS", hessian=TRUE, control=list(trace=1)))
fit3$par
##
system.time(fit4 <-
            .Call("test_cox_tvc3",
                  list(time=lung3$time,event=lung3$event,x=lung3$x,beta=beta),package="rstpm2"))
fit4
system.time(fit <- coxph(Surv(time, event) ~ x + tt(x), data=lung3,
             tt=function(x,t,...) x*log(t)))
(coxph(Surv(time, event) ~ x, data=lung3))

fit0 <- pstpm2(Surv(time, event) ~ x, data=lung3, smooth.formula=~s(log(time)))
fit <- pstpm2(Surv(time, event) ~ 1, data=lung3, smooth.formula=~s(log(time))+s(log(time),by=x),sp.init=c(1,1))
fit2 <- pstpm2(Surv(time, event) ~ x, data=lung3, smooth.formula=~s(log(time))+x:log(time),sp.init=1)
summary(fit2)
##
i <- -(1:10)
coef1 <- coef(fit)[i]
statistic <- as.numeric(coef1 %*% solve(vcov(fit)[i,i]) %*% coef1)
1-pchisq(statistic,length(i))
##
fit0 <- stpm2(Surv(time, event) ~ x, data=lung3, df=3)
fit <- stpm2(Surv(time, event) ~ 1, data=lung3, tvc.formula=~x:ns(log(time),df=3))
fit2 <- stpm2(Surv(time, event) ~ x, data=lung3, tvc.formula=~x:log(time), df=4)
summary(fit2)
anova(fit2,fit0)
##
## Wald test by hand
i <- 5:7
coef1 <- coef(fit)[i]
statistic <- as.numeric(coef1 %*% solve(vcov(fit)[i,i]) %*% coef1)
1-pchisq(statistic,length(i))
##
system.time(fit <- pstpm2(Surv(time, event) ~ x, data=lung3, smooth.formula=~s(log(time))+log(time):x,sp.init=1))
src <- '
    List largs = as<List>(args);
    vec time   = as<vec>(largs["time"]); // length n
    vec event  = as<vec>(largs["event"]); // length n
    mat X      = as<mat>(largs["X"]); // design matrix, n*c
    vec beta   = as<vec>(largs["beta"]); // length c+1
    int k      = as<int>(largs["k"]); // column to use for tvc
    int n      = time.size();
    int c      = X.n_cols;
    double llike = 0.0, lsum;
    vec eta;
    vec eta0 = X * beta(span(0,c-1));
    for (int i=0; i<n; ++i) {
      if (event(i) == 1 && (i==0 || time(i) != time(i-1))) {
	eta = eta0(span(i,n-1)) + beta(c)*log(time(i))*X(span(i,n-1),k);
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
    mat X      = as<mat>(largs["X"]); // one covariate, length n
    vec beta   = as<vec>(largs["beta"]); // length 2
    int k      = as<int>(largs["k"]); // column to use for tvc
    int n      = time.size();
    int c      = X.n_cols;
    vec grad(beta.size(),fill::zeros);
    vec lsum(beta.size(),fill::zeros);
    vec eta0 = X * beta(span(0,c-1));
    vec risk;
    mat Xrisk;
    for (int i=0; i<n; ++i) {
      if (event(i) == 1 && (i==0 || time(i) != time(i-1))) {
	risk = exp(eta0(span(i,n-1)) + beta(c)*log(time(i))*X(span(i,n-1),k));
        Xrisk = X(span(i,n-1), span::all);
        Xrisk.each_col() %= risk;
        lsum(span(0,c-1)) = sum(Xrisk, 0)/sum(risk);
        lsum(c) = sum(log(time(i))*Xrisk(span::all,k))/sum(risk);
        for (int j=i; j<n && time(j)==time(i) && event(j)==1; ++j) {
	  grad(span(0,c-1)) += X(j,span::all) - lsum(span(0,c-1));
	  grad(c) += X(j,k)*log(time(i)) - lsum(c);
        }
      }
    }
    return wrap(grad);
'
gr <- cxxfunction(signature(args="list"), src, plugin='RcppArmadillo', rcpp_inc)

.Call("test_cox_tvc2",with(lung3,list(time=time,event=event,x=x,beta=beta)),package="rstpm2")
fn0(with(lung3,list(time=time,event=event,x=x,beta=beta)))
fn(with(lung3,list(time=time,event=event,X=cbind(x),beta=beta,k=0)))
.Call("test_cox_tvc2_grad",with(lung3,list(time=time,event=event,x=x,beta=beta)),package="rstpm2")
gr0(with(lung3,list(time=time,event=event,x=x,beta=beta)))
gr(with(lung3,list(time=time,event=event,X=cbind(x),beta=beta,k=0)))

system.time(fit4 <-
            .Call("test_cox_tvc3",
                  list(time=lung3$time,event=lung3$event,x=lung3$x,beta=beta),package="rstpm2"))
fit4$coef
sqrt(diag(solve(fit4$hessian)))

system.time(fit2 <- optim(par=beta,
                          fn = function(beta) 
                              -fn0(list(time=lung3$time,event=lung3$event,x=lung3$x,beta=beta)),
                          gr = function(beta) 
                              -gr0(list(time=lung3$time,event=lung3$event,x=lung3$x,beta=beta)),
                          method="BFGS", hessian=TRUE, control=list(trace=1)))
fit2
system.time(fit3 <- optim(par=beta,
                          fn = function(beta) 
                              -fn(list(time=lung3$time,event=lung3$event,X=cbind(lung3$x),beta=beta,k=0)),
                          gr = function(beta) 
                              -gr(list(time=lung3$time,event=lung3$event,X=cbind(lung3$x),beta=beta,k=0)),
                          method="BFGS", hessian=TRUE, control=list(trace=1)))
fit3$par
sqrt(diag(solve(fit3$hessian)))
## // Stata output:
## ------------------------------------------------------------------------------
##           _t |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
## -------------+----------------------------------------------------------------
## main         |
##            x |   1.920236   .9687835     1.98   0.047     .0214556    3.819017
## -------------+----------------------------------------------------------------
## tvc          |
##            x |  -.3795621   .1682714    -2.26   0.024     -.709368   -.0497561
## ------------------------------------------------------------------------------


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
