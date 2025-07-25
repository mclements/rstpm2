## package.skeleton(name="rstpm2", path="c:/usr/src/R", force=TRUE, namespace=TRUE, code_files="pm2-3.R")
## Local Windows setup:
## Rtools.bat
## R CMD INSTALL --html "c:/usr/src/R/rstpm2/pkg"
## R CMD build "c:/usr/src/R/rstpm2/pkg"
## R CMD INSTALL --build "c:/usr/src/R/rstpm2/pkg"
## R CMD CHECK "c:/usr/src/R/rstpm2/pkg"
##
## Local Ubuntu setup:
## R CMD INSTALL --html ~/src/R/rstpm2/pkg --library=~/R/x86_64-pc-linux-gnu-library/2.12
## R CMD build ~/src/R/rstpm2/pkg
## R CMD build --binary ~/src/R/rstpm2/pkg
##
## testPackage <- TRUE
## if (testPackage) {
##   require(splines)
##   require(survival)
##   require(bbmle)
## }

## bug report
library(rstpm2)
fit <- aft(Surv(rectime,censrec==1)~hormon,data=brcancer,df=3) # ok
fit <- aft(Surv(rectime,censrec==1)~1,data=brcancer,df=3) # now ok

fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,df=3) # ok

## Issue with S3/S4 conflicts
library(rstpm2)
fit <- aft(Surv(rectime,censrec==1)~hormon,data=brcancer,df=3) # ok
inherits(fit, "mle2")
class(fit)
methods("predictnl")
plot(fit, type="surv", newdata=data.frame(hormon=1))

## bug report from Alex
library(rstpm2)
fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,df=3)
update(fit, formula. = . ~ . + x2)

## nsx for cure
library(rstpm2)
FUyear = seq(0,20,len=301)[-1]
n1 = nsx(log(FUyear), Boundary.knots = log(c(5,15)), knots=log(c(10,14)), cure=TRUE)
matplot(FUyear,n1,type="l")

## bug in tinyplot
library(tinyplot)
x=y=0:1
conf.low=y-0.5
conf.high=y+0.5
par(mfrow=1:2)
plt(y~x,ymin=conf.low,ymax=conf.high,type="ribbon")
plt(y~x,ymin=conf.low,ymax=conf.high,type="ribbon",xlim=c(0,0.5),ylim=c(0,0.5))

plt(y~x,type="l") |> names()
plt(y~x,type="l",xlim=c(0,0.5),ylim=c(0,0.5)) # ok
plt(y~x,ymin=conf.low,ymax=conf.high,type="l")
plt(y~x,ymin=conf.low,ymax=conf.high,type="l",xlim=c(0,0.5),ylim=c(0,0.5))


## Bug report
library(rstpm2)
colon3 <- data.frame(
    sex = ifelse(colon$sex == "Female", 0, 1),
    status = ifelse(colon$status %in% c("Dead: cancer", "Dead: other"), 1, ifelse(colon$status == "Alive", 0, NA)),
    persontime = as.numeric(colon$exit - colon$dx))
modvc <- stpm2(Surv(persontime, status) ~ sex, df = 3, tvc = list(sex = 3),
               data = colon3)

plot(modvc,
     newdata = data.frame(sex = 0),
     type = "hr",
     var = "sex",
     ci = TRUE, rug = FALSE,
     las = 1, ylim = c(0.8,1.5), lwd = 2,
     main = "Hazard Ratio (TVC)",
     xaxt = "n",
     xlab = "Time (years)")
axis(1, at = c(seq(0, 365*20, 365)), labels = c(0:20))
abline(v=365,lty=2)
df = data.frame(sex=0,persontime=365*c(1,2,5))
predict(modvc, type="hr", var="sex", newdata=df, se.fit=TRUE,full=TRUE) |>
    transform(labels=sprintf("%.3f\n(%.3f, %.3f)", Estimate, lower, upper)) |>
    with({abline(v=persontime,lty=2)
        use.upper=c(TRUE,FALSE,TRUE)
        text(x=persontime,upper*use.upper+lower*(!use.upper),labels=labels,
             pos=ifelse(use.upper,3,1))})


library(rstpm2)
library(broom)
test1 <- stpm2(Surv(rectime, censrec == 1) ~ hormon, data = brcancer, df = 3)
tidy(test1)
tidy(test1, conf.int=TRUE, exponentiate=TRUE)


model.matrix(~factor(x)-1, data.frame(x=1:2))

library(rstpm2)
aft(Surv(rectime,censrec==1)~hormon,data=rstpm2::brcancer,df=3) # ok
aft(Surv(rectime,censrec==1)~factor(hormon),data=rstpm2::brcancer,df=3) # ok
aft(Surv(rectime,censrec==1)~hormon+factor(hormon),data=rstpm2::brcancer,df=3) # ok
aft(Surv(rectime,censrec==1)~factor(hormon),data=rstpm2::brcancer,df=3,tvc=list(hormon=2)) # ok (drops second spline term)
aft(Surv(rectime,censrec==1)~hormon,data=rstpm2::brcancer,df=3,tvc=list(hormon=2)) # ok (drops hormon)
aft(Surv(rectime,censrec==1)~1,data=rstpm2::brcancer,df=3,tvc=list(hormon=2)) # ok

## marginaleffects
library(rstpm2)
library(marginaleffects)
test1 <- stpm2(Surv(rectime, censrec == 1) ~ hormon * x3, data = brcancer, df = 3)
nd <- data.frame(rectime=1000, hormon=c(0,1), x3=50)
pred1 <- predict(m, type = "surv", newdata = nd)
pred2 <- predictions(m, type = "surv", newdata = nd)
all(pred1 == pred2$estimate)

pred3 <- predict(m, newdata=nd, type="meansurv")
pred4 <- avg_predictions(m, type="surv", newdata=nd, by="x3")
pred3 == pred4$estimate

## Bug fix: aft with factor
library(rstpm2)
summary(fit1 <- aft(Surv(rectime,censrec==1)~hormon,data=rstpm2::brcancer,df=4))
summary(fit2 <- aft(Surv(rectime,censrec==1)~factor(hormon),data=rstpm2::brcancer,df=4))
vcov(fit1) |> cov2cor() |> "rownames<-"(NULL) |> "colnames<-"(NULL)
vcov(fit2) |> cov2cor() |> "rownames<-"(NULL) |> "colnames<-"(NULL)

## Bug report: aft with left truncation
library(rstpm2)
brcancer$start <- 0
fit <- aft(Surv(start, rectime,censrec==1)~hormon,data=brcancer,df=4) # now okay:)
brcancer$start <- 1
fit <- aft(Surv(start, rectime,censrec==1)~hormon,data=brcancer,df=4) # now okay:)

length(rstpm2:::lhs(Surv(start, rectime, censrec == 1) ~ hormon))

## test for predict(..., type="lpmatrixD")
library(rstpm2)
fit = aft(Surv(rectime,censrec==1)~hormon,data=brcancer)
plot(fit, newdata=data.frame(hormon=1))
lines(fit, newdata=data.frame(hormon=0), lty=2)


## test for predict(..., type="lpmatrixD")
library(rstpm2)
fit = stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer)
m = predict(fit, newdata=data.frame(hormon=1), grid=TRUE, type="lpmatrix")
m2 = predict(fit, newdata=data.frame(hormon=1), grid=TRUE, type="lpmatrixD")

## check that offset works
library(rstpm2)
stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer)
stpm2(Surv(rectime,censrec==1)~hormon+offset(hormon),data=brcancer)

library(survival)
library(rstpm2)
brcancer = transform(rstpm2::brcancer, id=1:(nrow(brcancer)/2))
survreg(Surv(rectime,censrec==1)~hormon+cluster(id),data=brcancer)
survreg(Surv(rectime,censrec==1)~hormon,data=brcancer)


## Andreas's bug report
library(rstpm2)
two_states <- function(model, ...) {
  transmat = matrix(c(NA,1,NA,NA),2,2,byrow=TRUE)
  rownames(transmat) <- colnames(transmat) <- c("Initial","Final")
  rstpm2::markov_msm(list(model), ..., trans = transmat)
}
## ERROR
death = gsm(Surv(time,status)~factor(rx), 
            data=survival::colon, subset=(etype==2), df=3)
ts = two_states(death, newdata=data.frame(rx = levels(survival::colon$rx)), 
                t = seq(0,2500, length = 51))
sts <- standardise(ts)
as.data.frame(sts)
##
## OK
death = gsm(Surv(time,status)~factor(rx) + age, 
            data=survival::colon, subset=(etype==2), df=3)
ts = two_states(death, newdata=data.frame(rx="Obs",
                                          age = c(65, 70, 75)), 
                t = seq(0,2500, length = 51))
sts <- standardise(ts)
as.data.frame(sts)


## predictions
library(rstpm2)
fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer, df=3, tvc=list(hormon=2))
plot(fit, type="hr", newdata=data.frame(hormon=0), var="hormon")
predict(fit, newdata=data.frame(hormon=0, rectime=c(250,500,1000)), type="hr", var="hormon", full=TRUE, se.fit=TRUE)
##
fit0 <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer, df=3)
anova(fit0,fit)

Y = 10
T = 1000
poisson.test(Y,T)
## log(theta)~normal
exp(log(Y)+c(-1.96,1.96)/sqrt(Y))
## theta~gamma
qgamma(c(0.025,0.975), c(Y,Y+1), 1)
Ys = rgamma(1e6, Y+1, 1)/T
quantile(Ys, c(0.025,0.975))
rates = exp(rnorm(1e5, log(Y), 1/sqrt(Y)))/T
quantile(rates, c(0.025,0.975))

## left truncated
library(rstpm2)
brcancer2 <- transform(brcancer,
                       startTime=ifelse(hormon==0,rectime/2,0))
##debug(stpm2)
summary(fit <- stpm2(Surv(startTime,rectime,censrec==1)~hormon,data=brcancer2, df=3))
tmp=predict(fit, newdata=data.frame(hormon=0), type="hazard", grid=TRUE, full=TRUE, se.fit=TRUE)


library(rstpm2)
library(biostat3)
library(survival)
colon = transform(biostat3::colon,
                  male=0+(sex=="Male"),
                  Unknown=0+(stage=="Unknown"),
                  Localised=0+(stage=="Localised"),
                  Regional=0+(stage=="Regional"),
                  Distant=0+(stage=="Distant"))
localised = subset(colon, stage=="Localised")
system.time(fit1 <- aft(Surv(surv_mm, status=="Dead: cancer") ~ I(age-50) + male, df=4,
                        data=localised, mixture=TRUE, tvc.integrated=TRUE,
                        control=list(reltol=1e-12), tvc=list(male=2)))

system.time(fit2 <- aft(Surv(surv_mm, status=="Dead: cancer") ~ I(age-50) + male, df=4,
                        data=localised, tvc=list(male=2)))



library(Rcpp)
library(RcppArmadillo)
sourceCpp(file = "~/src/R/rstpm2/inst/vintegrateRcpp2.cpp")
vintegrate = function(f, lower, upper, ..., subdivisions = 100L, rel.tol=.Machine$double.eps^0.25,
                       abs.tol=rel.tol) {
    f <- match.fun(f)
    ff <- function(x) f(x, ...)
    ny = max(length(lower), length(upper))
    if (ny>1 && length(lower)==1) lower = rep(lower,ny)
    if (ny>1 && length(upper)==1) upper = rep(upper,ny)
    limit <- as.integer(subdivisions)
    if (limit < 1L || (abs.tol <= 0 && rel.tol < max(50 * .Machine$double.eps, 
                                                     5e-29))) 
        stop("invalid parameter values")
    if (any(is.na(lower) | is.na(upper)))
        stop("a limit is NA or NaN")
    if (!(all(is.finite(lower)) || all(!is.finite(lower))))
        stop("lower bounds are a mix of finite and infinite values")
    if (!(all(is.finite(upper)) || all(!is.finite(upper))))
        stop("upper bounds are a mix of finite and infinite values")
    if (all(is.finite(lower)) && all(is.finite(upper))) {
        res = vdqagsRcpp(ff, lower, upper, rel.tol, abs.tol, limit, ny)
        res$value = drop(res$value)
        res$abs.err = drop(res$abs.err)
        res$abserr = NULL
    } else {
        if (all(is.finite(lower))) {
            inf <- 1L
            bound <- lower
        }
        else if (all(is.finite(upper))) {
            inf <- -1L
            bound <- upper
        }
        else {
            inf <- 2L
            bound <- rep(0,ny) # ignored:)
        }
        res = vdqagiRcpp(ff, bound, inf, rel.tol, abs.tol, limit, ny)
        res$value = drop(res$value)
        res$abserr = drop(res$abserr)
    }
    res
}
        
vintegrate(exp, -Inf, 1)$value - exp(1) # ok (scalar)
vintegrate(function(x) exp(c(x[1],2*x[2])), -Inf, c(1,1))$value - c(exp(1),exp(2)/2)
vintegrate(exp, -Inf, 1:2)$value - exp(1:2) # ok (one function and two bounds)
vintegrate(function(x) exp(x*1:2), -Inf, 1:2)$value - c(exp(1),exp(4)/2) # ok (vector function with vector bounds)
vintegrate(function(x) exp(2*x), -Inf, 1:2)$value - exp(c(2,4))/2
vintegrate(function(x) c(exp(x[1]),log(x[2])), 0, c(1,1))$value - c(exp(1)-1, -1)
vintegrate(function(x) c(exp(x[1]),log(2*x[2])), c(0,1), 2:3)$value

integrate(exp, 0, 2)
integrate(function(x) log(2*x), 1, 3)

vintegrate(function(x) c(exp(x[1]),log(2*x[2])), c(0,1), 2:3)

test_vdqags()
test_vdqags()$result- c(exp(1)-1, (exp(2)-1)/2, -1)

tm1 <- system.time(for (i in 1:10000) integrate(function(x) exp(2*x), -Inf, 0))
tm2 <- system.time(tmp <- vintegrate(function(x) exp(2*x), -Inf, rep(0,10000), subdivisions = 10))
tm1/tm2
head(tmp$abs.err)

vintegrate(function(x) log(2*x), 0, 10)
unclass(integrate(function(x) log(2*x), 0, 10))
tm1 <- system.time(for (i in 1:10000) integrate(function(x) log(2*x), 0, 10))
tm2 <- system.time(tmp <- vintegrate(function(x) log(2*x), rep(0,10000), 10, subdivisions = 20))
tm1/tm2
lapply(tmp,head)

vrdqk21Rcpp(function(x) log(c(x,x)), 0, 1)
## Rcpp::List vrdqk15Rcpp(Rcpp::Function f, double boun, int inf, double a, double b);
           ## inf    - int
           ##          indicating the kind of integration range involved
           ##          inf = 1 corresponds to  (bound,+infinity),
           ##          inf = -1            to  (-infinity,bound),
           ##          inf = 2             to (-infinity,+infinity).
vrdqk15Rcpp(f=function(x) exp(c(x,x)), boun=1.0, inf=-1, a=0, b=1)
integrate(exp, -Inf, 1)

(tmp <- vdqagsRcpp(function(x) log(c(x,2*x)), 0, 1, 1e-6, 1e-6, 15, 2))

vdqagsRcpp(function(x) exp(c(x,2*x)), 0, 1, 1e-6, 1e-6, 50, 2)
c(exp(1)-1,(exp(2)-1)/2)

vdqagsRcpp(log, 0, 1, 1e-6, 1e-6, 50, 1)
vdqagsRcpp(function(x) log(rep(x,2)), 0, 1, 1e-6, 1e-6, 50, 2)


## tm1 <- system.time(vdqagsRcpp(function(x) log(rep(x,10000)), 0, 1, 1e-6, 1e-6, 50, 10000)) # FAILS
tm1 <- system.time(vdqagsRcpp(function(x) log(rep(x,10000)), 0, 1, 1e-6, 1e-6, 10, 10000))
tm2 <- system.time(for (i in 1:10000) integrate(log, 0, 1))
tm2/tm1 # Only 3 times faster -- and does not scale for memory:(

## Nicely behaved functions are faster:
tm1 <- system.time(vdqagsRcpp(function(x) exp(rep(x,10000)), 0, 1, 1e-6, 1e-6, 10, 10000))
tm2 <- system.time(for (i in 1:10000) integrate(exp, 0, 1))
tm2/tm1 # 30 times faster:)




## testing for voptimize
library(rstpm2)
voptimize(function(x) {y=dnorm(x, c(0,1), log=TRUE); cat("x:",x,"\n","y:",y,"\n"); y}, lower=c(-5,-5), upper=c(5,5), maximum=TRUE, tol=1e-5)

set.seed(12345)
z=rnorm(1e4)
tm1 <- system.time(tmp <- voptimize(function(x) dnorm(x, z, log=TRUE), lower=c(-10,-10), upper=c(10,10), maximum=TRUE, tol=1e-5))
tm2 <- system.time(for (i in 1:length(z)) optimize(function(x) dnorm(x, z[i], log=TRUE), lower=c(-10,-10), upper=c(10,10), maximum=TRUE, tol=1e-5))
tm2/tm1

## Conditional distributions
set.seed(12345)
grs_variance = 0.68 ## variance on the log scale or on the frailty scale??
other_variance = 1.14
grs_log_mean = -grs_variance/2
other_log_mean = -other_variance/2
grs_log_sd = sqrt(grs_variance)
other_log_sd = sqrt(grs_variance)
## conditional on the sum (on the log scale), what is the GRS?
f1g = function(g,k,mu1=0,sigma1=1,mu2=0,sigma2=1)
    integrate(function(y1) g(y1)*dnorm(y1,mu1,sigma1)*dnorm(k-y1,mu2,sigma2), -Inf, Inf)$value
f1c = function(y1,k,mu1=0,sigma1=1,mu2=0,sigma2=1) {
    denom = f1g(function(x) 1, k, mu1, sigma1, mu2, sigma2)
    dnorm(y1,mu1,sigma1)*dnorm(k-y1,mu2,sigma2)/denom
}
dbinorm_total = function(x,k,mu1=0,sigma1=1,mu2=0,sigma2=1)
    dnorm(x,mu1,sigma1)*dnorm(k-x,mu2,sigma2)
integrate(f1c, -Inf, Inf, k=2,mu1=1,sigma1=2) # check: value==1
print(muc<-integrate(function(x) x*f1c(x,k=2,mu1=1,sigma1=2), -Inf, Inf))
print(varc<-integrate(function(x) (x-muc$value)^2*f1c(x,k=2,mu1=1,sigma1=2), -Inf, Inf))
y1=seq(-10,10,length=1001)
plot(y1, f1c(y1, 2,mu1=1,sigma1=2), type="l")
lines(y1,dnorm(y1,muc$value,sqrt(varc$value)),col="blue")
abline(v=1.8,lty=2)
##
set.seed(12345)
r1 = rnorm(1e5, muc$value, sqrt(varc$value))
r2 = 2-r1
c(mean(r1),sd(r1))
c(mean(r2),sd(r2))
plot(density(r1),xlim=c(-10,10))
lines(density(r2),lty=2)

library(Rcpp)
code = "#include <R_ext/Applic.h>
#include <Rcpp.h>
struct sbinorm_total { double k, mu1, sd1, mu2, sd2, norm, mean; };
void dbinorm_total(double *x, int n, void *ex) {
  sbinorm_total *p = (sbinorm_total *) ex;
  for (int i=0; i<n; ++i) 
    x[i] = R::dnorm(x[i], p->mu1, p->sd1, 0)*R::dnorm(p->k-x[i], p->mu2, p->sd2, 0);
}
void dbinorm_total_expected(double *x, int n, void *ex) {
  sbinorm_total *p = (sbinorm_total *) ex;
  // for (int i=0; i<n; ++i) Rprintf(\"%g\\n\",x[i]);
  for (int i=0; i<n; ++i) 
    x[i] = x[i]*R::dnorm(x[i], p->mu1, p->sd1, 0)*R::dnorm(p->k-x[i], p->mu2, p->sd2, 0)/p->norm;
}
void dbinorm_total_var(double *x, int n, void *ex) {
  sbinorm_total *p = (sbinorm_total *) ex;
  // for (int i=0; i<n; ++i) Rprintf(\"%g\\n\",x[i]);
  for (int i=0; i<n; ++i) 
    x[i] = (x[i]-p->mean)*(x[i]-p->mean)*R::dnorm(x[i], p->mu1, p->sd1, 0)*R::dnorm(p->k-x[i], p->mu2, p->sd2, 0)/p->norm;
}
// [[Rcpp::export]]
Rcpp::NumericVector test(Rcpp::NumericVector x, double k=2.0, double mu1=0.0, double sd1=1.0, double mu2=0.0, double sd2=1.0) {
  double * xx = &x[0];
  int n = x.size();
  sbinorm_total p{k, mu1, sd1, mu2, sd2, 0.0, 0.0};
  dbinorm_total(&xx[0], n, (void *) &p);
  Rcpp::NumericVector out(xx, xx+n);
  return out;
}
// [[Rcpp::export]]
double binorm_total_norm(double k=2.0, double mu1=0.0, double sd1=1.0, double mu2=0.0, double sd2=1.0) {
  sbinorm_total p{k, mu1, sd1, mu2, sd2};
  int inf=2, neval=0, ier=0, limit=1000, last=0, lenw;
  double epsabs = 1.0e-6, epsrel = 1.0e-6, result = 0.0,
    abserr = 0.0, bound = 0.0;
  lenw = 4 * limit;
  int *iwork =   (int *) R_alloc(limit, sizeof(int));
  double *work = (double *) R_alloc(lenw,  sizeof(double));
  Rdqagi(dbinorm_total, (void *) &p, &bound, &inf,
            &epsabs, &epsrel,
            &result, &abserr, &neval, &ier,
            &limit, &lenw, &last,
            iwork, work);
  return result;
}
// [[Rcpp::export]]
double binorm_total_expected(double k=2.0, double mu1=0.0, double sd1=1.0, double mu2=0.0, double sd2=1.0) {
  int inf=2, neval=0, ier=0, limit=1000, last=0, lenw;
  double epsabs = 1.0e-4, epsrel = 1.0e-4, result = 0.0, 
    abserr = 0.0, bound = 0.0;
  lenw = 4 * limit;
  int *iwork =   (int *) R_alloc(limit, sizeof(int));
  double *work = (double *) R_alloc(lenw,  sizeof(double));
  double norm=binorm_total_norm(k,mu1,sd1,mu2,sd2);
  sbinorm_total p2{k, mu1, sd1, mu2, sd2, norm, 0.0};
  Rdqagi(dbinorm_total_expected, (void *) &p2, &bound, &inf,
            &epsabs, &epsrel,
            &result, &abserr, &neval, &ier,
            &limit, &lenw, &last,
            iwork, work);
  return result;
}
// [[Rcpp::export]]
Rcpp::NumericVector binorm_total_sd(double k=2.0, double mu1=0.0, double sd1=1.0, double mu2=0.0, double sd2=1.0) {
  int inf=2, neval=0, ier=0, limit=1000, last=0, lenw;
  double epsabs = 1.0e-6, epsrel = 1.0e-6, result = 0.0, 
    abserr = 0.0, bound = 0.0;
  lenw = 4 * limit;
  int *iwork =   (int *) R_alloc(limit, sizeof(int));
  double *work = (double *) R_alloc(lenw,  sizeof(double));
  sbinorm_total p{k, mu1, sd1, mu2, sd2, 0.0, 0.0};
  Rdqagi(dbinorm_total, (void *) &p, &bound, &inf,
            &epsabs, &epsrel,
            &result, &abserr, &neval, &ier,
            &limit, &lenw, &last,
            iwork, work);
  p.norm = result;
  Rdqagi(dbinorm_total_expected, (void *) &p, &bound, &inf,
            &epsabs, &epsrel,
            &result, &abserr, &neval, &ier,
            &limit, &lenw, &last,
            iwork, work);
  p.mean = result;
  Rdqagi(dbinorm_total_var, (void *) &p, &bound, &inf,
            &epsabs, &epsrel,
            &result, &abserr, &neval, &ier,
            &limit, &lenw, &last,
            iwork, work);
  Rcpp::NumericVector out{p.mean, std::sqrt(result)};
  return out;
}
// [[Rcpp::export]]
Rcpp::NumericMatrix test2(Rcpp::NumericVector k, double mu1=0.0, double sd1=1.0, double mu2=0.0, double sd2=1.0) {
  Rcpp::NumericMatrix out(k.size(),3);
  Rcpp::NumericVector y(2);
  for (int i=0; i<k.size(); ++i) {
    y = binorm_total_sd(k[i], mu1, sd1, mu2, sd2);
    out(i,0) = k[i];
    out(i,1) = y[0];
    out(i,2) = y[1];
  }
  return out;
}
"
sourceCpp(code=code)
test(0:1,2,1,2,0,1)
dbinorm_total(0:1, 2, 1, 2)
binorm_total_norm(2,1,2,0,1)
f1g(function(x) 1, 2, mu1=1, sigma1=2, mu2=0, sigma2=1)
binorm_total_expected(2,1,2)
print(mu <- integrate(function(x) x*f1c(x,k=2,mu1=1,sigma1=2), -Inf, Inf))
binorm_total_var(2,1,2)
integrate(function(x) (x-mu$value)^2*f1c(x,k=2,mu1=1,sigma1=2), -Inf, Inf)

binorm_total_sd(k=2,mu1=1,sd1=2,mu2=0,sd2=1)
binorm_total_sd(k=0,mu1=0,sd1=1,mu2=0,sd2=1)


print(mu <- integrate(function(x) x*f1c(x,k=0,mu1=1,sigma1=2), -Inf, Inf))
binorm_total_var(0,1,2)
integrate(function(x) (x-mu$value)^2*f1c(x,k=0,mu1=1,sigma1=2), -Inf, Inf)

system.time(test2(seq(-6,6,length=1e5), mu1=grs_log_mean, sd1=grs_log_sd, mu2=other_log_mean, sd2=other_log_sd))
system.time(tmp <- test2(seq(-6,6,0.1), mu1=grs_log_mean, sd1=grs_log_sd, mu2=other_log_mean, sd2=other_log_sd))
tmp

y = rnorm(1e5, grs_log_mean+other_log_mean, sqrt(grs_log_sd^2+other_log_sd^2))
range(y)

system.time(for (i in 1:1000) binorm_total_var(i))

## Integration
f = function(x) {cat(sort(x), "\n"); log(x)}
integrate(f, 0, 1)
f = (function() {n=0; xs = c(); function(x) {n <<- n+1; xs <<- c(xs,x); log(x)}})()
integrate(f, 0, 1)
get("n",environment(f))
length(get("xs",environment(f)))
plot(density(get("xs",environment(f)),from=0,to=1))
##
stats::integrate(log,0,1)
rmutil::int(function(x) log(x),0,1, eps=1e-10)+1
system.time(for (i in 1:100) integrate(log,0,1))

myintegrate = function(f, a, b, control=list()) {
    control = modifyList(list(tol=1e-8,maxn=20, debug=FALSE),
                         control)
    xi = c(.995657163025808080735527280689003,
            .973906528517171720077964012084452,
            .930157491355708226001207180059508,
            .865063366688984510732096688423493,
            .780817726586416897063717578345042,
            .679409568299024406234327365114874,
            .562757134668604683339000099272694,
            .433395394129247190799265943165784,
            .294392862701460198131126603103866,
            .14887433898163121088482600112972,
            0.)
    xi = c(-xi, xi[10:1])
    w = c(.066671344308688137593568809893332,
	    .149451349150580593145776339657697,
	    .219086362515982043995534934228163,
	    .269266719309996355091226921569469,
	    .295524224714752870173892994651338)
    w = c(w, w[5:1]) / 2
    wk = c(.011694638867371874278064396062192,
            .03255816230796472747881897245939,
            .05475589657435199603138130024458,
            .07503967481091995276704314091619,
            .093125454583697605535065465083366,
            .109387158802297641899210590325805,
            .123491976262065851077958109831074,
            .134709217311473325928054001771707,
            .142775938577060080797094273138717,
            .147739104901338491374841515972068,
            .149445554002916905664936468389821)
    wk = c(wk, wk[10:1]) / 2
    step = function(a,b,n) {
        if (control$debug)
            cat(sprintf("(%g,%g), n=%i\n",a,b,n))
        fx = f ( a + (b - a ) * ( xi + 1)/2)
        qb1 = (b - a ) * sum(w * fx[seq.int(2,20,2)])
        qa1 = (b - a ) * sum(wk * fx)
        xi2 = c( ( xi - 1)/2 , ( xi + 1)/2 )
        fx2 = f ( a + (b - a )*( xi2 +1)/2)
        qb2 = (b - a ) * sum(c( w , w ) * fx2[ c(seq(2,20,2) , seq(23,41,2))]) / 2
        qa2 = (b - a ) * sum(c( wk , wk ) * fx2) / 2
        int = qa2
        qab1 = qa1 - qb1 ; qaa = qa2 - qa1 ;
        qab2 = qa2 - qb2 ; qbb = qb2 - qb1 ;
        noise = 21 * .Machine$double.eps * (b - a ) * max ( abs( fx ))
        if (abs( qab1 ) <= noise) qab1 = 0
        if (abs( qab2 ) <= noise) qab2 = 0
        if (abs( qaa ) <= noise) qaa = 0
        if (abs( qbb ) <= noise) qbb = 0
        if ((qab1==0 || (0<=qab2/qab1 && qab2/qab1 <= 1)) &&
            (qaa==0 || abs(qaa) <= abs(qbb))) {
            serr = if(qaa==0) 0 else qab2*qaa/(qbb-qaa)
        } else serr = 2*tol
        attr(int,"serr") = serr
        attr(int,"err") = abs(serr)
        if (abs(serr)<=control$tol || n>=control$maxn) {
            if (n>=control$maxn) warning("Exceeded the maximum recursive depth")
            return(int)
        } else {
            m = a+(b-a)/2
            int1 = step(a, m, n+1)
            int2 = step(m, b, n+1)
            int = int1+int2
            attr(int,"serr") = attr(int1,"serr")+attr(int2,"serr")
            attr(int,"err") = attr(int1,"err")+attr(int2,"err")
            return(int)
        }
    }
    out = step(a,b,0)
    out+attr(out,"serr")
}
myintegrate(log,0,1,control=list(tol=1e-12)) + 1

## Based on gsl/integration/{qag.c,qk.c}: see integration.cpp and vintegration*.cpp
library(Rcpp)
library(RcppArmadillo)
sourceCpp(file="~/src/R/rstpm2/inst/vintegrationRcpp.cpp")
vintegrate = function(f, a, b, abstol=1e-8, reltol=1e-8, maxit=1000) {
    ff = function(x) f(a+(b-a)*x)
    res = vintegrationRcpp(ff, abstol, reltol, maxit)
    res$value = drop(res$value*(b-a))
    res$abserr = drop(res$abserr*(b-a))
    res
}
vintegrate(log, 0, 1)
vintegrate(function(x) log(c(x,x,x,x)), 0, 1)
system.time(vintegrate(function(x) rep(log(x),100), 0, 1))
system.time(for (i in 1:100) integrate(log, 0, 1))

sourceCpp(file="~/src/R/rstpm2/inst/integrateRcpp.cpp")
integrateRcpp = function(f, lower, upper, ..., rel.tol=.Machine$double.eps^0.25, abs.tol=rel.tol, subdivisions = 100L) {
    f <- match.fun(f)
    ff <- function(x) f(x, ...)
    limit <- as.integer(subdivisions)
    if (limit < 1L || (abs.tol <= 0 && rel.tol < max(50 * .Machine$double.eps, 
                                                     5e-29))) 
        stop("invalid parameter values")
    stopifnot(length(lower) == 1, length(upper) == 1)
    dqagsRcpp(f, lower, upper, abs.tol, rel.tol, subdivisions)
}
integrateRcpp(log, 0, 1)
rdqk21Rcpp(function(x) {cat(x,"\n"); exp(x)}, 0, 1)
integrate(function(x) {cat(x,"\n"); exp(x)}, 0, 1)

sourceCpp(file="~/src/R/rstpm2/inst/vintegrateRcpp.cpp")
vintegrateRcpp = function(f, lower, upper, ..., rel.tol=.Machine$double.eps^0.25, abs.tol=rel.tol, subdivisions = 100L) {
    f <- match.fun(f)
    ff <- function(x) f(lower+(upper-lower)*x, ...)
    fmid = ff(0.5)
    limit <- as.integer(subdivisions)
    if (limit < 1L || (abs.tol <= 0 && rel.tol < max(50 * .Machine$double.eps, 
                                                     5e-29))) 
        stop("invalid parameter values")
    stopifnot(length(lower) == 1, length(upper) == 1)
    res = vdqagsRcpp(ff, lower, upper, abs.tol, rel.tol, subdivisions)
    
}
vintegrateRcpp(log, 0, 1)
integrate(function(x) {cat(x,"\n"); exp(x)}, 0, 1)

rdqk21Rcpp(function(x) {cat(x,"\n"); log(x)}, 0, 1)
vrdqk21Rcpp(function(x) {cat(x,"\n"); log(x)}, 0, 1)
vrdqk21Rcpp(function(x) log(c(x,x)), 0, 1)


## Based on https://github.com/jacobwilliams/quadpack/blob/master/src/quadpack_generic.F90
dqk21 = function(f, a, b, ...) {
    ## n = 11
    xgk = c(.995657163025808080735527280689003,
            .973906528517171720077964012084452,
            .930157491355708226001207180059508,
            .865063366688984510732096688423493,
            .780817726586416897063717578345042,
            .679409568299024406234327365114874,
            .562757134668604683339000099272694,
            .433395394129247190799265943165784,
            .294392862701460198131126603103866,
            .14887433898163121088482600112972,
            0.)
    ## xgk = c(-xgk, xgk[10:1])
    wgk = c(.011694638867371874278064396062192,
            .03255816230796472747881897245939,
            .05475589657435199603138130024458,
            .07503967481091995276704314091619,
            .093125454583697605535065465083366,
            .109387158802297641899210590325805,
            .123491976262065851077958109831074,
            .134709217311473325928054001771707,
            .142775938577060080797094273138717,
            .147739104901338491374841515972068,
            .149445554002916905664936468389821)
    ## wgk = c(wgk, wgk[10:1]) / 2
    wg = c(.066671344308688137593568809893332,
	    .149451349150580593145776339657697,
	    .219086362515982043995534934228163,
	    .269266719309996355091226921569469,
           .295524224714752870173892994651338)
    ## wg = c(wg, wg[5:1]) / 2
    fv1 <- fv2 <- vector("numeric",10)
    epmach = .Machine$double.eps
    centr = 0.5*(a+b)
    hlgth = 0.5*(b-a)
    dhlgth = abs(hlgth)
    ##
    resg = 0.0
    fc = f(centr, ...)
    resk = wgk[11]*fc
    Resabs = abs(resk)
    for (j in 1:5) {
        jtw = 2L*j
        absc = hlgth*xgk[jtw]
        fval1 = f(centr - absc, ...)
        fval2 = f(centr + absc, ...)
        fv1[jtw] = fval1
        fv2[jtw] = fval2
        fsum = fval1 + fval2
        resg = resg + wg[j]*fsum
        resk = resk + wgk[jtw]*fsum
        Resabs = Resabs + wgk[jtw]*(abs(fval1) + abs(fval2))
    }
    for (j in 1:5) {
        jtwm1 = 2L*j - 1L
        absc = hlgth*xgk[jtwm1]
        fval1 = f(centr - absc, ...)
        fval2 = f(centr + absc, ...)
        fv1[jtwm1] = fval1
        fv2[jtwm1] = fval2
        fsum = fval1 + fval2
        resk = resk + wgk[jtwm1]*fsum
        Resabs = Resabs + wgk[jtwm1]*(abs(fval1) + abs(fval2))
    }
    reskh = resk*0.5
    Resasc = wgk[11]*abs(fc - reskh)
    for (j in 1:10) 
        Resasc = Resasc + wgk[j]*(abs(fv1[j] - reskh) + abs(fv2[j] - reskh))
    Result = resk*hlgth
    Resabs = Resabs*dhlgth
    Resasc = Resasc*dhlgth
    Abserr = abs((resk - resg)*hlgth)
    if (Resasc != 0.0 && Abserr != 0.0) 
            Abserr = Resasc*min(1.0, (200.0*Abserr/Resasc)^1.5)
    if (Resabs > uflow/(50.0*epmach))
        Abserr = max((epmach*50.0)*Resabs, Abserr)
    c(Result=Result, Abserr=Abserr, Resabs=Resabs, Resasc=Resasc)
}
dqage = function(f, a, b, Epsabs=1e-10, Epsrel=1e-10, Limit=15L, ...) {
    epmach = .Machine$double.eps
    uflow = .Machine$double.xmin
    Key = 2
    Alist <- Blist <- Elist <- Rlist <- vector("numeric", Limit)
    Iord <- vector("integer", Limit)
    Ier = 0L
    Neval = 0L
    Last = 0L
    Result = 0.0
    Abserr = 0.0
    Alist[1] = a
    Blist[1] = b
    Rlist[1] = 0.0
    Elist[1] = 0.0
    Iord[1] = 0L
    if (Epsabs <= 0.0 && Epsrel < max(50.0*epmach, 0.5e-28)) Ier = 6L
    if (Ier != 6) then
    Neval = 0L
    res = dqk21(f, a, b, ...)
    Last = 1L
    Rlist[1] = res$Result
    Elist[1] = res$Abserr
    Iord[1] = 1L
    errbnd = max(Epsabs, Epsrel*abs(Result))
    if (Abserr <= 50.0*epmach*defabs && Abserr > errbnd) Ier = 2L
    if (Limit == 1) Ier = 1L
    if (! (Ier != 0 || (Abserr <= errbnd && Abserr != resabs) 
        || Abserr == 0.0)) {
        errmax = Abserr
        maxerr = 1L
        area = Result
        errsum = Abserr
        nrmax = 1L
        iroff1 = 0L
        iroff2 = 0L
        for (Last in 2:Limit) {
            a1 = Alist[maxerr]
            b1 = 0.5*(Alist[maxerr] + Blist[maxerr])
            a2 = b1
            b2 = Blist[maxerr]
            ## call dqk21(f, a1, b1, area1, error1, resabs, defab1)
            ## call dqk21(f, a2, b2, area2, error2, resabs, defab2)
            ## c(Result=Result, Abserr=Abserr, Resabs=Resabs, Resasc=Resasc)
            res1 = dqk21(f, a1, b1, ...)
            area1 = res1$Result
            error1 = res1$Abserr
            defab1 = res1$Resasc
            res2 = dqk21(f, a2, b2, ...)
            area2 = res2$Result
            error2 = res2$Abserr
            defab2 = res2$Resasc
            Neval = Neval + 1L
            area12 = area1 + area2
            erro12 = error1 + error2
            errsum = errsum + erro12 - errmax
            area = area + area12 - Rlist[maxerr]
            if (defab1 != error1 && defab2 != error2) {
                if (abs(Rlist[maxerr] - area12) <= 0.1e-4*abs(area12) 
                            && erro12 >= 0.99*errmax) iroff1 = iroff1 + 1L
                if (Last > 10 && erro12 > errmax) iroff2 = iroff2 + 1L
            }
            Rlist[maxerr] = area1
            Rlist[Last] = area2
            errbnd = max(Epsabs, Epsrel*abs(area))
            if (errsum > errbnd) {
                if (iroff1 >= 6L || iroff2 >= 20L) Ier = 2L
                if (Last == Limit) Ier = 1L
                if (max(abs(a1), abs(b2)) <=
                    (1.0 + 100.0*epmach)*(abs(a2) + 1000.0*uflow)) Ier = 3L
            }
            if (error2 > error1) {
                Alist[maxerr] = a2
                Alist[Last] = a1
                Blist[Last] = b1
                Rlist[maxerr] = area2
                Rlist[Last] = area1
                Elist[maxerr] = error2
                Elist[Last] = error1
            } else {
                Alist[Last] = a2
                Blist[maxerr] = b1
                Blist[Last] = b2
                Elist[maxerr] = error1
                Elist[Last] = error2
            }
        }
    }
}

myintegrate = function(f, a, b, control=list()) {
    control = modifyList(list(tol=1e-8,maxn=20, debug=FALSE),
                         control)
    step = function(a,b,n) {
        if (control$debug)
            cat(sprintf("(%g,%g), n=%i\n",a,b,n))
        fx = f ( a + (b - a ) * ( xi + 1)/2)
        qb1 = (b - a ) * sum(w * fx[seq.int(2,20,2)])
        qa1 = (b - a ) * sum(wk * fx)
        xi2 = c( ( xi - 1)/2 , ( xi + 1)/2 )
        fx2 = f ( a + (b - a )*( xi2 +1)/2)
        qb2 = (b - a ) * sum(c( w , w ) * fx2[ c(seq(2,20,2) , seq(23,41,2))]) / 2
        qa2 = (b - a ) * sum(c( wk , wk ) * fx2) / 2
        int = qa2
        qab1 = qa1 - qb1 ; qaa = qa2 - qa1 ;
        qab2 = qa2 - qb2 ; qbb = qb2 - qb1 ;
        noise = 21 * .Machine$double.eps * (b - a ) * max ( abs( fx ))
        if (abs( qab1 ) <= noise) qab1 = 0
        if (abs( qab2 ) <= noise) qab2 = 0
        if (abs( qaa ) <= noise) qaa = 0
        if (abs( qbb ) <= noise) qbb = 0
        if ((qab1==0 || (0<=qab2/qab1 && qab2/qab1 <= 1)) &&
            (qaa==0 || abs(qaa) <= abs(qbb))) {
            serr = if(qaa==0) 0 else qab2*qaa/(qbb-qaa)
        } else serr = 2*tol
        attr(int,"serr") = serr
        attr(int,"err") = abs(serr)
        if (abs(serr)<=control$tol || n>=control$maxn) {
            if (n>=control$maxn) warning("Exceeded the maximum recursive depth")
            return(int)
        } else {
            m = a+(b-a)/2
            int1 = step(a, m, n+1)
            int2 = step(m, b, n+1)
            int = int1+int2
            attr(int,"serr") = attr(int1,"serr")+attr(int2,"serr")
            attr(int,"err") = attr(int1,"err")+attr(int2,"err")
            return(int)
        }
    }
    out = step(a,b,0)
    out+attr(out,"serr")
}
myintegrate(log,0,1,control=list(tol=1e-12)) + 1


mintegrate = function(f, a, b, control=list(), ...) {
    stopifnot(length(f(a))==length(a),
              length(a) == length(b))
    control = modifyList(list(tol=1e-8,maxn=20, debug=FALSE),
                         control)
    xi = c(.995657163025808080735527280689003,
            .973906528517171720077964012084452,
            .930157491355708226001207180059508,
            .865063366688984510732096688423493,
            .780817726586416897063717578345042,
            .679409568299024406234327365114874,
            .562757134668604683339000099272694,
            .433395394129247190799265943165784,
            .294392862701460198131126603103866,
            .14887433898163121088482600112972,
            0.)
    xi = c(-xi, xi[10:1])
    w = c(.066671344308688137593568809893332,
	    .149451349150580593145776339657697,
	    .219086362515982043995534934228163,
	    .269266719309996355091226921569469,
	    .295524224714752870173892994651338)
    w = c(w, w[5:1]) / 2
    wk = c(.011694638867371874278064396062192,
            .03255816230796472747881897245939,
            .05475589657435199603138130024458,
            .07503967481091995276704314091619,
            .093125454583697605535065465083366,
            .109387158802297641899210590325805,
            .123491976262065851077958109831074,
            .134709217311473325928054001771707,
            .142775938577060080797094273138717,
            .147739104901338491374841515972068,
            .149445554002916905664936468389821)
    wk = c(wk, wk[10:1]) / 2
    width = b-a
    olda = a
    oldb = b
    ff = function(x) f(olda+x*width, ...)
    ## fmid = ff(0.5)
    step = function(a,b,n,side) {
        if (control$debug)
            cat(sprintf("(%g,%g), n=%i\n",a,b,n))
        fx = t(sapply(xi, function(xii) ff ( a + (b - a ) * ( xii + 1)/2)))
        qb1 = (b - a ) * colSums(w * fx[seq.int(2,20,2),])
        qa1 = (b - a ) * colSums(wk * fx)
        xi2 = c( ( xi - 1)/2 , ( xi + 1)/2 )
        fx2 = t(sapply(xi2, function(xii2) ff ( a + (b - a )*( xii2 +1)/2)))
        qb2 = (b - a ) * colSums(c( w , w ) * fx2[ c(seq(2,20,2) , seq(23,41,2)),]) / 2
        qa2 = (b - a ) * colSums(c( wk , wk ) * fx2) / 2
        int = qa2
        qab1 = qa1 - qb1 ; qaa = qa2 - qa1 ;
        qab2 = qa2 - qb2 ; qbb = qb2 - qb1 ;
        noise = 21 * .Machine$double.eps * (b - a ) * max ( abs( fx ))
        qab1 = ifelse(abs( qab1 ) <= noise, 0, qab1)
        qab2 = ifelse(abs( qab2 ) <= noise, 0, qab2)
        qaa = ifelse(abs( qaa ) <= noise, 0, qaa)
        qbb = ifelse(abs( qbb ) <= noise, 0, qbb)
        serr = ifelse((qab1==0 | (0<=qab2/qab1 & qab2/qab1 <= 1)) &
                      (qaa==0 | abs(qaa) <= abs(qbb)),
               ifelse(qaa==0, 0, qab2*qaa/(qbb-qaa)),
               rep(2*tol,length(qaa)))
        attr(int,"serr") = serr
        attr(int,"err") = abs(serr)
        if (all(abs(serr)<=control$tol) || n>=control$maxn) {
            if (n>=control$maxn && side==TRUE) warning("Exceeded the maximum recursive depth")
            return(int)
        } else {
            m = a+(b-a)/2
            int1 = step(a, m, n+1, TRUE)
            int2 = step(m, b, n+1, FALSE)
            int = int1+int2
            attr(int,"serr") = attr(int1,"serr")+attr(int2,"serr")
            attr(int,"err") = attr(int1,"err")+attr(int2,"err")
            return(int)
        }
    }
    out = step(0,1,0, TRUE)
    out+attr(out,"serr")
}
system.time(mintegrate(function(x) c(rep(exp(x),10000), log(x)),c(rep(0,10000),1),c(rep(1,10000),2),control=list(tol=1e-12)) - c(rep(exp(1)-1,10000), 2*log(2) - 1))
system.time(mintegrate(function(x) c(rep(exp(x),10), log(x)),c(rep(0,10),1),c(rep(1,10),2),control=list(tol=1e-12)) - c(rep(exp(1)-1,10), 2*log(2) - 1))
system.time({for (i in 1:10000) integrate(log,0,1); integrate(exp,0,1)})

system.time(mintegrate(function(x) c(exp(x[1:10000]), log(x[11])),c(rep(0,10000),1),c(rep(1,10000),2),control=list(tol=1e-12)) - c(rep(exp(1)-1,10000), 2*log(2) - 1))
mintegrate(function(x) c(exp(x[1:10]), log(x[11])),c(rep(0,10),1),c(rep(1,10),2),control=list(tol=1e-12)) - c(rep(exp(1)-1,10), 2*log(2) - 1)

myintegrate(function(x) rep(1,length(x)),0,1) - 1
myintegrate(function(x) x,0,1) - 1/2
myintegrate(function(x) x^2,0,1) - 1/3
myintegrate(function(x) x^3,0,1) - 1/4
myintegrate(function(x) x*log(x),0,1,control=list(tol=1e-6)) + 1/4
myintegrate(function(x) 1/x,0,1) - 1
myintegrate(function(x) 1/sqrt(x),0,1) - 2
integrate(function(x) 1/sqrt(x),0,1)
integrate(function(x) 1/sqrt(x),0,1)$value -2


## AFT model fitting with linear constraints
library(rstpm2)
library(biostat3)
library(survival)
colon = transform(biostat3::colon,
                  male=0+(sex=="Male"),
                  Unknown=0+(stage=="Unknown"),
                  Localised=0+(stage=="Localised"),
                  Regional=0+(stage=="Regional"),
                  Distant=0+(stage=="Distant"))
localised = subset(colon, stage=="Localised")

system.time(fit1 <- stpm2(Surv(surv_mm, status=="Dead: cancer") ~ I(age-50) + factor(stage), df=4, data=colon))
pred=predict(fit1,newdata=transform(colon,stage="Localised"),type='meansurvdiff',grid=TRUE,full=TRUE,se.fit=TRUE,exposed=function(data) transform(data, stage="Distant"))

## Get the penalty matrix for the baseline splines
system.time(fit1 <- aft(Surv(surv_mm, status=="Dead: cancer") ~ I(age-50) + male, df=4, data=localised, mixture=TRUE, tvc.integrated=TRUE, control=list(reltol=1e-12), tvc=list(male=2))) # 10.3s
plot(fit1, newdata=data.frame(age=50, male=1))
system.time(fit1a <- aft(Surv(surv_mm, status=="Dead: cancer") ~ I(age-50), df=4, data=localised, mixture=TRUE, tvc.integrated=FALSE, control=list(constrOptim=TRUE, reltol=1e-12), tvc=list(male=2))) # 13.5s
fit1
fit1@args$negll(coef(fit1))
fit1@args$negll(-coef(fit1))
fit1@args$negll(coef(fit1a))
fit1@args$gradient(coef(fit1))
fit1@args$gradient(coef(fit1a))
fit1@args$gradient(-coef(fit1))

ui = fit1@args$ui
negll = \(coef) fit1@args$negll(coef,add.penalties=FALSE)
grad = \(coef) fit1@args$gradient(coef, add.penalties=FALSE)
## negll = fit1@args$negll
## grad = fit1@args$gradient
system.time(fit1b <- constrOptim(coef(fit1)*1.05, f=negll,
                                 grad=grad, ui=ui, ci=0, mu=1e-10, outer.eps=1e-10,
                                 control=list(reltol=1e-12)))
fit1b
print(fit1b$value,digits=12)
fit1@args$gradient(fit1b$par)
##
library(alabama)
system.time(al1 <- alabama::auglag(coef(fit1)*1.05, fn=negll, gr=grad,
                                   hin=function(coef) ui %*% coef,
                                   hin.jac=function(coef) ui,
                                   control.optim=list(rel.tol=1e-12),
                                   control.outer=list(eps=1e-10,
                                                      method="nlminb")))
al1
##
library(nloptr)
system.time(nl1 <- nloptr::auglag(coef(fit1)*1.05, fn=negll,
                                  hin=function(coef) ui %*% coef,
                                  gr=grad,
                                  hinjac=function(coef) ui,
                                  localsolver = "LBFGS",
                                  localtol=1e-10))
nl1
fit1@args$gradient(nl1$par)
##
par(mfrow=c(2,2))
plot(coef(fit1)/al1$par)
plot(coef(fit1a)/al1$par)
plot(fit1b$par/al1$par)
plot(nl1$par/al1$par)

aft(Surv(surv_mm, status=="Dead: cancer") ~ I(age-50), df=8, data=localised, mixture=TRUE, tvc=list(male=2), reltol=1e-12, init=-coef(fit1)) # wrong:(

print(fit1@args$negll(coef(fit1)), digits=12)
print(fit1@args$negll(coef(fit1a)), digits=12)
print(fit1@args$negll(fit1b$par), digits=12)
print(fit1@args$negll(coef(fit1b)), digits=12)

coef(fit1)
coef(fit1a)

aft(Surv(surv_mm, status=="Dead: cancer") ~ I(age-50) + male, df=2, data=localised)
aft_mixture(Surv(surv_mm, status=="Dead: cancer") ~ I(age-50) + male, df=2, data=localised)
aft_integrated(Surv(surv_mm, status=="Dead: cancer") ~ I(age-50) + male, df=2, data=localised)
survreg(Surv(surv_mm, status=="Dead: cancer") ~ I(age-50) + male, data=localised)

aft(Surv(surv_mm, status=="Dead: cancer") ~ I(age-50) + male, df=2, data=localised)
aft(Surv(surv_mm, status=="Dead: cancer") ~ I(age-50) + male + bhazard(bhaz), df=2, data=transform(localised,bhaz=1e-4))

aft_mixture(Surv(surv_mm, status=="Dead: cancer") ~ I(age-50) + male, df=2, data=localised, mixture=TRUE, cure.formula=~1)
aft(Surv(surv_mm, status=="Dead: cancer") ~ I(age-50) + male, df=2, data=localised, mixture=TRUE)
aft_integrated(Surv(surv_mm, status=="Dead: cancer") ~ I(age-50) + male, df=2, data=localised, mixture=TRUE, cure.formula=~1)
aft(Surv(surv_mm, status=="Dead: cancer") ~ I(age-50) + male, df=2, data=localised, mixture=TRUE, tvc.integrated=TRUE)

fit1 = aft_mixture(Surv(surv_mm, status=="Dead: cancer") ~ I(age-50) + male, df=8, data=localised, mixture=TRUE, tvc=list(male=2))
fit2 = aft(Surv(surv_mm, status=="Dead: cancer") ~ I(age-50) + male, df=8, data=localised, mixture=TRUE, tvc=list(male=2))
fit3 = aft(Surv(surv_mm, status=="Dead: cancer") ~ I(age-50) + male, df=8, data=localised, mixture=TRUE, tvc.integrated=TRUE, tvc=list(male=2))


aft(Surv(surv_mm, status=="Dead: cancer") ~ I(age-50) + male, df=2, data=localised) # ok
aft(Surv(surv_mm, status=="Dead: cancer") ~ I(age-50) + male + bhazard(0*age), df=2, data=localised) # ok
aft(Surv(surv_mm, status=="Dead: cancer") ~ I(age-50) + male + bhazard(0*age), df=8, data=localised, mixture=TRUE, tvc=list(male=2)) # fails


fit2 = aft_integrated(Surv(surv_mm, status=="Dead: cancer") ~ I(age-50) + male, df=2, data=localised)
fit3 = survreg(Surv(surv_mm, status=="Dead: cancer") ~ I(age-50) + male, data=localised)


## Examples for the second AFT paper
library(rstpm2)
library(biostat3)
library(survival)
colon = transform(biostat3::colon,
                  male=0+(sex=="Male"),
                  Unknown=0+(stage=="Unknown"),
                  Localised=0+(stage=="Localised"),
                  Regional=0+(stage=="Regional"),
                  Distant=0+(stage=="Distant"))
localised = subset(colon, stage=="Localised")

## Story
## Some sanity checks -> move to tests
## Introduce the dataset: simulated data based on a population-based cancer registry
## We could also use the brcancer dataset
## Time-varying acceleration factors: stage or sex?

## What about left truncation? We could use population-based data for prostate cancer incidence or time to next PSA test or biopsy. Both of these analyses could use the Stockholm PSA and Biopsy Register.
## What about recurrent events? This would require either a marginal model, a copula model or a random effects model. This could use Stan or ADMB or Benjamin's library. For the last option, we could simplify Benjamin's VAJointSurv library -- and then think about other extensions.
## What about relative survival / excess hazard modelling? This would be a modest extension.

## show that aft(..., df=2) and survreg(...) give the same results for the beta coefficients
fit0 = aft(Surv(surv_mm, status=="Dead: cancer") ~ I(age-50) + male, df=2, data=localised)
fit1 = aft_mixture(Surv(surv_mm, status=="Dead: cancer") ~ I(age-50) + male, df=2, data=localised)
fit2 = aft_integrated(Surv(surv_mm, status=="Dead: cancer") ~ I(age-50) + male, df=2, data=localised)
fit3 = survreg(Surv(surv_mm, status=="Dead: cancer") ~ I(age-50) + male, data=localised)
fit0
fit1
fit2
fit3

## The following models also give the same estimates...
fit0 = aft(Surv(surv_mm, status=="Dead: cancer") ~ I(age-50) + male, data=localised)
fit1 = aft_mixture(Surv(surv_mm, status=="Dead: cancer") ~ I(age-50) + male, data=localised)
fit2 = aft_integrated(Surv(surv_mm, status=="Dead: cancer") ~ I(age-50) + male, data=localised)
fit0
fit1
fit2
par(mfrow=c(3,2))
plot(fit0, type="surv", newdata=data.frame(age=50, male=0), ylim=c(0.7,1), main="Females") # ok
plot(fit0, type="surv", newdata=data.frame(age=50, male=1), ylim=c(0.7,1), main="Males")   # ok
plot(fit1, type="surv", newdata=data.frame(age=50, male=0), ylim=c(0.7,1), main="Females") # ok
plot(fit1, type="surv", newdata=data.frame(age=50, male=1), ylim=c(0.7,1), main="Males")   # ok
plot(fit2, type="surv", newdata=data.frame(age=50, male=0), ylim=c(0.7,1), main="Females") # ok
plot(fit2, type="surv", newdata=data.frame(age=50, male=1), ylim=c(0.7,1), main="Males")   # ok
## plot(fit0, type="surv", newdata=data.frame(age=50, male=0:1)) ## issue with plotting multiple rows

## plot survival
library(survival)
summary(aft(Surv(surv_mm, status=="Dead: cancer") ~ male, data=localised))

## Kaplan-Meier estimators versus a fitted AFT model
par(mfrow=1:2)
summary(fit1 <- aft(Surv(surv_mm, status=="Dead: cancer") ~ male, df=4, data=localised))
plot(fit1, newdata=data.frame(male=1), type="surv",
     xlab="Time since cancer diagnosis (months)",
     main="Males",
     ylim=c(0.4,1))
lines(survfit(Surv(surv_mm, status=="Dead: cancer") ~ 1, data=localised,
              subset=(male==1)))
plot(fit1, newdata=data.frame(male=0), type="surv",
     xlab="Time since cancer diagnosis (months)",
     main="Females",
     ylim=c(0.4,1))
lines(survfit(Surv(surv_mm, status=="Dead: cancer") ~ 1, data=localised,
              subset=(male==0)))

## Female patients are older -> potential confounding
t.test(age~male, localised) 
plot(density(subset(localised,male==1)$age))
lines(density(subset(localised,male==0)$age), lty=2)

summary(fit0 <- aft(Surv(surv_mm, status=="Dead: cancer") ~ ns(age,df=2)+male, data=localised, df=4,
                    tvc=list(male=3), tvc.intercept=FALSE))
summary(fit1 <- aft_mixture(Surv(surv_mm, status=="Dead: cancer") ~ ns(age,df=2)+male, data=localised, df=4,
                    tvc=list(male=3), tvc.intercept=FALSE))
summary(fit2 <- aft_integrated(Surv(surv_mm, status=="Dead: cancer") ~ ns(age,df=2)+male, data=localised, df=4, tvc=list(male=3), tvc.intercept=FALSE))
par(mfrow=c(1,3))
plot(fit0, type="accfac", newdata=data.frame(age=70, male=0), ylim=c(0,2), var="male", main="aft")
plot(fit1, type="accfac", newdata=data.frame(age=70, male=0), ylim=c(0,2), var="male",
     main="aft_mixture")
plot(fit2, type="accfac", newdata=data.frame(age=70, male=0), ylim=c(0,2), var="male",
     main="aft_integrated")

## 
summary(fit0 <- aft(Surv(surv_mm, status=="Dead: cancer") ~ ns(age,df=2)+male+Distant, data=colon, df=4,
                    tvc=list(Distant=3), tvc.intercept=FALSE))
summary(fit1 <- aft_mixture(Surv(surv_mm, status=="Dead: cancer") ~ ns(age,df=2)+male+Distant, data=colon, df=4,
                    tvc=list(Distant=3), tvc.intercept=FALSE))
summary(fit2 <- aft_integrated(Surv(surv_mm, status=="Dead: cancer") ~ ns(age,df=2)+male+Distant, data=colon, df=4,
                    tvc=list(Distant=3), tvc.intercept=FALSE))
par(mfrow=c(1,3))
plot(fit0, type="accfac", newdata=data.frame(age=70, male=0, Distant=0), var="Distant", main="aft")
plot(fit1, type="accfac", newdata=data.frame(age=70, male=0, Distant=0), var="Distant",
     main="aft_mixture")
plot(fit2, type="accfac", newdata=data.frame(age=70, male=0, Distant=0), var="Distant",
     main="aft_integrated")

## 
summary(fit0 <- aft(Surv(surv_mm, status=="Dead: cancer") ~ ns(age,df=2)+male+Localised, data=colon, df=4,
                    tvc=list(Localised=3), tvc.intercept=FALSE))
summary(fit1 <- aft_mixture(Surv(surv_mm, status=="Dead: cancer") ~ ns(age,df=2)+male+Localised, data=colon, df=4,
                    tvc=list(Localised=3), tvc.intercept=FALSE))
summary(fit2 <- aft_integrated(Surv(surv_mm, status=="Dead: cancer") ~ ns(age,df=2)+male+Localised, data=colon, df=4,
                    tvc=list(Localised=3), tvc.intercept=FALSE))
par(mfrow=c(1,3))
plot(fit0, type="accfac", newdata=data.frame(age=70, male=0, Localised=0), var="Localised", main="aft")
plot(fit1, type="accfac", newdata=data.frame(age=70, male=0, Localised=0), var="Localised",
     main="aft_mixture")
plot(fit2, type="accfac", newdata=data.frame(age=70, male=0, Localised=0), var="Localised",
     main="aft_integrated")

## 
summary(fit0 <- aft(Surv(surv_mm, status=="Dead: cancer") ~ ns(age,df=2)+male+Distant, data=colon, df=4,
                    tvc=list(Distant=3, male=2), tvc.intercept=FALSE))
summary(fit1 <- aft_mixture(Surv(surv_mm, status=="Dead: cancer") ~ ns(age,df=2)+male+Distant, data=colon, df=4,
                    tvc=list(Distant=3, male=2), tvc.intercept=FALSE))
summary(fit2 <- aft_integrated(Surv(surv_mm, status=="Dead: cancer") ~ ns(age,df=2)+male+Distant, data=colon, df=4,
                    tvc=list(Distant=3,male=2), tvc.intercept=FALSE))
par(mfrow=c(1,3))
plot(fit0, type="accfac", newdata=data.frame(age=70, male=0, Distant=0), var="Distant", main="aft")
plot(fit1, type="accfac", newdata=data.frame(age=70, male=1, Distant=0), var="Distant",
     main="aft_mixture")
plot(fit2, type="accfac", newdata=data.frame(age=70, male=0, Distant=0), var="Distant",
     main="aft_integrated")

## Predictions for a tvc for one covariate does not change with a change in other covariates.
## I had expected that the aft() and aft_mixture formulations would lead to small differences in predictions. This is a pleasant property:).
predict(fit0, type="accfac", newdata=data.frame(age=70, male=0, Distant=0, surv_mm=50), var="Distant")
predict(fit0, type="accfac", newdata=data.frame(age=50, male=1, Distant=0, surv_mm=50), var="Distant")
predict(fit1, type="accfac", newdata=data.frame(age=70, male=0, Distant=0, surv_mm=50), var="Distant")
predict(fit1, type="accfac", newdata=data.frame(age=50, male=1, Distant=0, surv_mm=50), var="Distant")
predict(fit2, type="accfac", newdata=data.frame(age=70, male=0, Distant=0, surv_mm=50), var="Distant")
predict(fit2, type="accfac", newdata=data.frame(age=50, male=1, Distant=0, surv_mm=50), var="Distant")

## 
summary(fit0 <- aft(Surv(surv_mm, status=="Dead: cancer") ~ ns(age,df=2)+male+Distant, data=colon, df=4,
                    tvc=list(Distant=3, male=2), tvc.intercept=FALSE))
summary(fit1 <- aft_mixture(Surv(surv_mm, status=="Dead: cancer") ~ ns(age,df=2)+male+Distant, data=colon, df=4,
                    tvc=list(Distant=3, male=2), tvc.intercept=FALSE))
summary(fit2 <- aft_integrated(Surv(surv_mm, status=="Dead: cancer") ~ ns(age,df=2)+male+Distant, data=colon, df=4,
                    tvc=list(Distant=3,male=2), tvc.intercept=FALSE))
par(mfrow=c(1,3))
plot(fit0, type="sdiff", newdata=data.frame(age=70, male=0, Distant=0), var="Distant", main="aft")
plot(fit1, type="sdiff", newdata=data.frame(age=70, male=1, Distant=0), var="Distant",
     main="aft_mixture")
plot(fit2, type="sdiff", newdata=data.frame(age=70, male=0, Distant=0), var="Distant",
     main="aft_integrated")


summary(fit1 <- aft_mixture(Surv(surv_mm, status=="Dead: cancer") ~ male, data=localised, df=4, tvc=list(male=2), tvc.intercept=FALSE))
summary(fit2 <- aft_integrated(Surv(surv_mm, status=="Dead: cancer") ~ male, data=localised, df=4, tvc=list(male=2), tvc.intercept=FALSE))
par(mfrow=c(1,2))
plot(fit1, type="accfac", newdata=data.frame(age=70, male=0), ylim=c(0,2),
     exposed=function(data) transform(data,male=1)) # ok
plot(fit2, type="accfac", newdata=data.frame(age=70, male=0), ylim=c(0,2),
     exposed=function(data) transform(data,male=1)) # Does this make sense?


plot(fit1, newdata=data.frame(age=70, male=0), type="sdiff",
     ## exposed=function(data) transform(data, male=1),
     var="male",
     xlab="Time since cancer diagnosis (months)")

plot(fit2, newdata=data.frame(age=70, male=0), type="hr",
     ## exposed=function(data) transform(data, male=1),
     var="male",
     xlab="Time since cancer diagnosis (months)")



fit0 = aft(Surv(surv_mm, status=="Dead: cancer") ~ I(age-50), data=localised,
           tvc=list(male=3))
fit1 = aft_mixture(Surv(surv_mm, status=="Dead: cancer") ~ I(age-50), data=localised,
                   tvc=list(male=3))
fit2 = aft_integrated(Surv(surv_mm, status=="Dead: cancer") ~ I(age-50), data=localised,
                      tvc=list(male=3))
fit0
fit1
fit2
par(mfrow=c(1,3))
plot(fit0, type="accfac", newdata=data.frame(age=50, male=0), ylim=c(0,2),
     exposed=function(data) transform(data,male=1)) # ok
plot(fit1, type="accfac", newdata=data.frame(age=50, male=0), ylim=c(0,2),
     exposed=function(data) transform(data,male=1)) # ok
plot(fit2, type="accfac", newdata=data.frame(age=50, male=0), ylim=c(0,2),
     exposed=function(data) transform(data,male=1)) # ok

fit0 = aft(Surv(surv_mm, status=="Dead: cancer") ~ I(age-50) + male, data=localised,
           tvc.intercept=FALSE, tvc=list(male=2))
fit1 = aft_mixture(Surv(surv_mm, status=="Dead: cancer") ~ I(age-50) + male, data=localised,
                   tvc.intercept=FALSE, tvc=list(male=2))
fit2 = aft_integrated(Surv(surv_mm, status=="Dead: cancer") ~ I(age-50) + male, data=localised,
                      tvc.intercept=FALSE, tvc=list(male=2))
fit0
fit1
fit2
par(mfrow=c(1,3))
plot(fit0, type="accfac", newdata=data.frame(age=50, male=0), ylim=c(0,2),
     exposed=function(data) transform(data,male=1)) # ok
plot(fit1, type="accfac", newdata=data.frame(age=50, male=0), ylim=c(0,2),
     exposed=function(data) transform(data,male=1)) # ok
plot(fit2, type="accfac", newdata=data.frame(age=50, male=0), ylim=c(0,2),
     exposed=function(data) transform(data,male=1)) # ok

fit1 = aft_mixture(Surv(surv_mm, status=="Dead: cancer") ~ ns(age-50,df=3) + male, data=localised,
                   tvc=list(male=2))
fit2 = aft_integrated(Surv(surv_mm, status=="Dead: cancer") ~ ns(age-50,df=3) + male, data=localised,
                      tvc=list(male=2))
fit1
fit2
par(mfrow=c(1,2))
plot(fit1, type="accfac", newdata=data.frame(age=50, male=0,Unknown=0,Regional=0,Distant=0), ylim=c(0,2),
     exposed=function(data) transform(data,male=1)) # ok
plot(fit2, type="accfac", newdata=data.frame(age=50, male=0,Unknown=0,Regional=0,Distant=0), ylim=c(0,2),
     exposed=function(data) transform(data,male=1)) # ok

fit1 = aft_mixture(Surv(surv_mm, status=="Dead: cancer") ~ ns(age,df=3) + male + Unknown + Regional + Distant, data=colon, tvc=list(Distant=2))
fit2 = aft_integrated(Surv(surv_mm, status=="Dead: cancer") ~ ns(age,df=3) + male + Unknown + Regional + Distant, data=colon, tvc=list(Distant=2))
summary(fit1)
summary(fit2)
par(mfrow=c(1,2))
plot(fit1, type="accfac", newdata=data.frame(age=50, male=0,Unknown=0,Regional=0,Distant=0), ylim=c(0,2),
     exposed=function(data) transform(data,Distant=1)) # ok
plot(fit2, type="accfac", newdata=data.frame(age=50, male=0,Unknown=0,Regional=0,Distant=0), ylim=c(0,2),
     exposed=function(data) transform(data,Distant=1)) # ok



## Email from Breanna M
library(rstpm2)
## choose knots based on user-defined quantiles for the event times
aknots=quantile(log(subset(brcancer,censrec==1)$rectime),c(0.05,0.5,0.95))
nknots=length(aknots)
## fit using those knots
stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,
      smooth.formula=~ns(log(rectime),Boundary.knots=range(aknots),knots=aknots[-c(1,nknots)]))
## or (for shorter output)
myns = function(rectime) ns(log(rectime),Boundary.knots=range(aknots),knots=aknots[-c(1,nknots)])
stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,
      smooth.formula=~myns(rectime))



## Bug: gradients were different for aft and aft_mixture with tvc :(
fdiff = \(f,x,eps=1e-5)
    sapply(1:length(x),
           \(i) (f("[<-"(x,i,x[i]+eps)) - f("[<-"(x,i,x[i]-eps)))/2/eps)
fdiff(fit0@minuslogl, coef(fit0))
fit0@args$gradient(coef(fit0))
fdiff(fit1@minuslogl, coef(fit1))
fit1@args$gradient(coef(fit1))
##
fdiff(fit0@minuslogl, -coef(fit0))
fit0@args$gradient(-coef(fit0))
fdiff(fit1@minuslogl, -coef(fit1))
fit1@args$gradient(-coef(fit1))

summary(fit0) # ok - same fit as aft()
par(mfrow=1:2)
## debugonce(rstpm2:::plot.aft.base)
plot(fit0, type="surv", newdata=data.frame(age=50, male=0), main="Females") # ok
plot(fit0, type="surv", newdata=data.frame(age=50, male=1), main="Males")   # ok
## plot(fit0, type="surv", newdata=data.frame(age=50, male=0:1)) ## issue with plotting multiple rows

fit1 = aft(Surv(surv_mm, status=="Dead: cancer") ~ I(age-50) + male, data=localised, tvc = list(male=2))
fit1 = aft_mixture(Surv(surv_mm, status=="Dead: cancer") ~ I(age-50) + male, data=localised, tvc = list(male=2))
fit1 = aft_integrated(Surv(surv_mm, status=="Dead: cancer") ~ I(age-50) + male, data=localised, tvc = list(male=2))
summary(fit1)

par(mfrow=1:2)
plot(fit1, type="surv", newdata=data.frame(age=50, male=0), main="Females") # ok
plot(fit1, type="surv", newdata=data.frame(age=50, male=1), main="Males") # wrong-o!

plot(fit1, type="accfac", newdata=data.frame(age=50, male=0), var="male", ylim=c(0,3))


plot(fit1, type="accfac", newdata=data.frame(age=50, male=0),
     exposed=function(data) transform(data, male=1), ylim=c(0,3))

## Email from Paul S
library(rstpm2)
fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,df=3,tvc=list(hormon=2))
## Plot marginal HR
## Straight from plot function
plot(fit, type = "meanhr", newdata = transform(brcancer, hormon = 0), exposed = function(data) transform(data, hormon = 1))

## Using ggplot (need to manually extract predictions and plot)
library(ggplot2)
## debugonce(rstpm2:::predict.stpm2.base)
system.time(pred <- predict(fit,
                            newdata = transform(brcancer, hormon = 0),
                            grid = TRUE, full = TRUE, se.fit = TRUE, type = "meanhr", recent=TRUE,
                            exposed = function(data) transform(data, hormon = 1)))
system.time(pred2 <- predict(fit,
                             newdata = transform(brcancer, hormon = 0),
                             grid = TRUE, full = TRUE, se.fit = TRUE, type = "meanhr",
                             exposed = function(data) transform(data, hormon = 1)))
## Plot
ggplot(pred, aes(x = rectime, y = Estimate)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3) +
  ylim(0,2) +
  xlab('Time (days)') + ylab('HR') +
  theme_bw(base_size = 15)


## Email for episig
## install.packages("extRemes")
library(extRemes)
library(survival)
set.seed(12345)
n<-500000
x1<-rnorm(n)
x2<-rnorm(n)
error<-revd(n, loc = 0, scale = 1)
b0<- 1
b1<-0.5
b2<--1.2
c<-1
time<-exp(b0+b1*x1+b2*x2-c*error)
status<-rep(1,n)
survreg(Surv(time, status==1) ~ x1+x2,dist="exponential")


library(rstpm2)
object <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,df=3,tvc=list(hormon=2))
newdata = data.frame(hormon=0)
simulate(object, newdata=newdata) # currently returns just the times -- this will later return a data-frame
simulate(object, newdata=newdata, start=1500) # left-truncation
plot(object, newdata=newdata)

grep_call = function(name,x) {
    local_function = function(x)
        if(length(x)==1) x==name else any(sapply(x, local_function))
    which(sapply(x, local_function))
}
gsm_design = function(object, newdata, inflate=100) {
    stopifnot(inherits(object, "stpm2"),
              is.list(newdata),
              is.numeric(inflate),
              length(inflate) == 1)
    ## Assumed patterns:
    ## timeEffect := (ns|nsx)(log(timeVar),knots,Boundary.knots,centre=FALSE,derivs=(c(2,2)|c(2,1)))
    ## effect := timeEffect | otherEffect:timeEffect | timeEffect:otherEffect
    terms = attr(object@model.frame, "terms")
    factors = attr(terms, "factors")[-1,,drop=FALSE]
    variables = attr(terms, "variables")[-(1:2)]
    predvars = attr(terms, "predvars")[-(1:2)]
    indices = grep_call(object@timeVar, variables)
    if(length(indices)==0) stop("No timeVar in the formula -- unexpected error")
    index_time_variables = grep_call(object@timeVar, variables) # time variables
    index_time_effects = grep(object@timeVar,colnames(factors)) # components in the rhs with time variables
    ## I want to know how wide is each term
    nms = names(coef(object))
    term.labels = attr(terms, "term.labels")
    coef_index <-
        sapply(strsplit(nms, ":"), function(c) {
            if (length(c)>2)
                stop("current implementation only allows for main effects and two-way interactions")
            pmatchp = function(x, table)
                !is.na(pmatch(x, table))
            if(length(c)==1) {
                for (i in seq_along(term.labels)) {
                    t = term.labels[i]
                    ## browser()
                    if (pmatchp(t,c))
                        return(i)
                }
                if (c == "(Intercept)")
                    return(0)
                return(-1)
            }
            ## => length(c) == 2
            term.labels.split = strsplit(term.labels, ":")
            for (i in seq_along(term.labels)) {
                t = term.labels.split[[i]]
                if (all(pmatchp(t,c)))
                    return(i)
            }
            return(-1)
        })
    ## how to map from term.labels to factors to predvars?
    ## Can I replace each term with its predvars?
    parse_ns = function(mycall,x,index_time_effect) {
        df = length(c(mycall$knots, mycall$Boundary.knots)) - 1
        stopifnot(mycall[[1]] == quote(nsx) || mycall[[1]] == quote(ns),
                  length(mycall[[2]])>1,
                  mycall[[2]][[1]] == quote(log), # assumes log
                  mycall[[2]][[2]] == object@timeVar, # what about a scalar product or divisor?
                  is.null(mycall$deriv) || (mycall$derivs[1] == 2 && mycall$derivs[2] %in% 1:2),
                  mycall$centre == FALSE) # doesn't allow for centering
        cure = !is.null(mycall$derives) && all(mycall$derivs == c(2,1))
        time = object@args$time
        q_const = attr(nsx(log(mean(time)), knots=mycall$knots,
                           Boundary.knots=mycall$Boundary.knots,
                           intercept=mycall$intercept),
                       "q.const")
        ## browser()
        list(call = mycall,
             knots=mycall$knots,
             Boundary_knots=mycall$Boundary.knots,
             intercept=as.integer(mycall$intercept),
             gamma=coef(object)[which(coef_index %in% index_time_effect)],
             q_const = q_const,
             cure = as.integer(cure),
             x=x)
    }
    time = object@args$time
    newdata[[object@timeVar]] = mean(time) # NB: time not used
    Xp = predict(object, newdata=newdata, type="lpmatrix")
    index2 = which(coef_index %in% index_time_effects)
    etap = drop(Xp[, index2, drop=FALSE] %*% coef(object)[index2])
    list(type="gsm",
         link_name=object@args$link,
         tmin = min(time), # not currently used?
         tmax = max(time),
         inflate=as.double(inflate),
         etap=etap,
         log_time=TRUE,
         terms =
             lapply(index_time_effects,
                    function(i) {
                        j = which(factors[,i] != 0)
                        if (length(j)==1)
                            return(parse_ns(predvars[[j]], rep(1, nrow(newdata)), i))
                        else {
                            if(length(j)>3)
                                stop("Current implementation only allows for two-way interaction terms")
                            if (j[1] %in% index_time_variables)
                                return(parse_ns(predvars[[j[1]]], eval(predvars[[j[2]]], newdata), i))
                            else return(parse_ns(predvars[[j[2]]], eval(predvars[[j[1]]], newdata), i))
                        }
                    })
         )
}

library(rstpm2)
library(microsimulation)
object <- gsm(Surv(rectime,censrec==1)~hormon,data=brcancer,df=3,tvc=list(hormon=2))
newdata = data.frame(hormon=0:1)
inflate = 100
(design <- gsm_design(object,newdata))
replicate(10,.Call("test_read_gsm", design, PACKAGE="microsimulation"))


set.seed(1002)
simulate(fit1, nsim=10, newdata=nd)
simulate(fit1, newdata=data.frame(hormon=0:1))
set.seed(1002)
replicate(10,.Call("test_read_gsm", design, PACKAGE="rstpm2"))
set.seed(1002)
range(replicate(10,.Call("test_read_gsm", gsm_design(fit1,nd), PACKAGE="rstpm2")) -
      simulate(fit1, nsim=10, seed=1002, newdata=nd))
library(microsimulation)
set.seed(1002)
replicate(10,.Call("test_read_gsm", gsm_design(fit1,nd), PACKAGE="microsimulation"))

simulate.stpm2 <- function(object, nsim=1, seed=NULL,
                           newdata=NULL, lower=1e-6, upper=1e5, ...) {
    stopifnot(!is.null(newdata))
    if (!is.null(seed)) set.seed(seed)
    ## assumes nsim replicates per row in newdata
    n = nsim * nrow(newdata)
    newdata = newdata[rep(1:nrow(newdata), each=nsim), , drop=FALSE]
    e <- rexp(n)
    objective <- function(time) {
        newdata[[object@timeVar]] <- time
        predict(object, type="cumhaz", newdata=newdata) - e
    }
    vuniroot(objective, lower=rep(lower,length=n), upper=rep(upper,length=n))$root
}
setMethod("simulate", "stpm2", simulate.stpm2)
setMethod("simulate", "pstpm2", simulate.stpm2)



## Can we estimate the hessian externally?
library(rstpm2)
library(numDeriv)
fit = stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,df=4)
summary(fit)
fit@vcov = solve(numDeriv::hessian(fit@minuslogl, coef(fit)))
summary(fit)
## Can we also use the gradient?
(solve(numDeriv::jacobian(fit@args$gradnegll, coef(fit))) - solve(numDeriv::hessian(fit@minuslogl, coef(fit)))) |> range()
(solve(numDeriv::jacobian(fit@args$gradnegll, coef(fit))) - vcov(fit)) |> range()

## Can we improve on accuracy using nsxD? This would also include predictions using C++ code.

##
library(rstpm2)
m <- rstpm2::stpm2(Surv(time, status) ~ sex, data = lung)
predict(m, newdata = data.frame(sex = c(1,1), time = c(100,200)), type = "rmst", se.fit=TRUE) # ok
plot(m, newdata = data.frame(sex = 1), type="rmst")
plot(m, newdata = data.frame(sex = 1), var="sex", type="rmstdiff")


## Predictions for differences
## Can we return the predictions, gradients and covariance matrix and use those for differences?
## preddiff = pred1 - pred0
## graddiff = grad1 - grad0
## vardiff = graddiff^T Sigma graddiff
## Note that we can silently return attributes
d = structure(data.frame(a=1), attr1=1:100, class=c("Test","data.frame"))
print(d) # print.Test() does not exist, so this uses print.data.frame()
## What about gradients on transformed scales,
## with a transformation g and an inverse transformation G?
## var(eta1) = grad_eta1^T Sigma grad_eta1
## and pred1 = G(eta1)
## => var(pred1) = var(G(eta1)) = grad(G(eta1))^T Sigma grad(G(eta1))
## where grad1 = grad(G(eta1)) = G'(eta1) grad_eta1?? (chain rule)
## preddiff = pred1 - pred0
## graddiff = grad1 - grad0
## vardiff = graddiff^T Sigma graddiff



## add cure for the AFT models
library(rstpm2)
fit0 = aft(Surv(rectime,censrec==1)~hormon,data=brcancer,df=4)
fit = aft(Surv(rectime,censrec==1)~hormon,data=brcancer,df=4,cure=TRUE)
par(mfrow=1:2)
plot(fit0,newdata=data.frame(hormon=0),main="Without cure")
plot(fit,newdata=data.frame(hormon=0),main="With cure")
cov2cor(vcov(fit0))

## cure models paper
library(cuRe)
colonDC <- subset(cuRe::colonDC, stage %in% c("Regional","Distant"))
colonDC <- transform(colonDC, stage=factor(stage),
                     stageDistant=as.numeric(stage == "Distant"),
                     bhaz = general.haz(time = "FU",
                                        rmap = list(age="agedays", sex = "sex", year = "dx"),
                                        data = colonDC, ratetable = survexp.dk))
fit.lat <- stpm2(Surv(FUyear, status) ~ stageDistant + bhazard(bhaz),
                 data = colonDC, df = 4, cure = TRUE)
fit.lat.timevar <- stpm2(Surv(FUyear, status) ~ stageDistant + bhazard(bhaz),
                         data = colonDC, df = 4, cure = TRUE,
                         tvc.formula = ~nsx(log(FUyear), df=4, cure=TRUE):stageDistant)
fit.lat.timevar2 <- stpm2(Surv(FUyear, status) ~ stageDistant + bhazard(bhaz),
                         data = colonDC, df = 4, cure = TRUE, tvc = list(stageDistant = 4))

predict(fit.lat, newdata = data.frame(FUyear = 0, stageDistant = 0), type = "probcure", se.fit = TRUE)
predict(fit.lat, newdata = data.frame(FUyear = 2, stageDistant = 0), type = "uncured", se.fit = TRUE)
predict(fit.lat.timevar, newdata = data.frame(FUyear = 2, stageDistant = 0), type = "probcure", se.fit = TRUE)

plot(fit.lat, newdata = data.frame(stageDistant = 0), type = "probcure") 
plot(fit.lat, newdata = data.frame(stageDistant = 0), type = "surv") 
plot(fit.lat, newdata = data.frame(stageDistant = 0), type = "uncured") 

## solve A x = b
set.seed(12345)
A = matrix(rnorm(18),3,byrow=TRUE)
b = c(0,0,1)
x = qr.solve(A,b) # under-determined system
drop(A %*% x - b) # this is a solution -- noice!
qmat = qr.Q(qr(t(A)), complete=TRUE)[, -(1:3), drop=FALSE]
A %*% qmat %*% rnorm(3)  # should be close to zero...
A %*% qmat

library(rstpm2)
library(splines)
print(ns1 <- nsx(1:10,df=4,intercept=TRUE))
## check how to calculate q.const
Bstar <- splineDesign(sort(c(rep(1,4),attr(ns1,"knots"),rep(10,4))),
                      x=c(1,10),derivs=c(2,2))
qr.Q(qr(t(Bstar)), complete=TRUE)[,-(1:2), drop=FALSE] - attr(ns1, "q.const") # OK
## Now, solve: Bstar x = c(0,0,1)
Bstar <- splineDesign(sort(c(rep(1,4),attr(ns1,"knots"),rep(10,4))),
                      x=c(1,10,10),derivs=c(2,2,1))
b = c(0,0,1)
x = qr.solve(Bstar,b) # magic sauce
drop(Bstar %*% x) - b # OK
qmat = qr.Q(qr(t(Bstar)), complete=TRUE)[,-(1:3), drop=FALSE]
theta = rnorm(3)
drop(Bstar %*% (drop(qmat %*% theta) + x)) - b # OK
Bstar %*% qmat %*% theta + Bstar %*% x - b # OK
## transformation: theta -> B %*% (drop(qmat %*% theta) + x) = B %*% qmat %*% theta + B %*% x
## nsx(..., derivs, target=rep(0,length(derivs))) for target=b
## returns an attribute for offset=x
## For fixed design matrix, we get N = B %*% qmat and C = B %*% x, so that theta -> N %*% theta + C
## How will this affect the linear predictor?
## Rather than N(t)^T theta, we will have N(t)^T theta + C
## For the derivative, we get N' = B' %*% qmat and C' = B' %*% x,
## so that we have N'(t)^T theta + C'
xs = seq(1,10,l=301)
B <- splineDesign(sort(c(rep(1,4),attr(ns1,"knots"),rep(10,4))),
                      x=xs)
N = B %*% qmat
C = B %*% x
theta = rnorm(3)
plot(xs, N %*% theta + C, type="l")
## Plot of derivatives...
B <- splineDesign(sort(c(rep(1,4),attr(ns1,"knots"),rep(10,4))),
                      x=xs, derivs=1)
N = B %*% qmat
C = B %*% x
theta = rnorm(3)
plot(xs, N %*% theta + C, type="l")
## How to predict outside of the boundaries?

## test for smoothpwc
library(rstpm2)
library(dplyr)
relative_survival <- function(model1, smoothpwc, ...) {
    transmat = matrix(c(NA,1,2,
                        NA,NA,NA,
                        NA,NA,NA),3,3,byrow=TRUE)
    rownames(transmat) <- colnames(transmat) <- c("Initial","Cause-specific death","Other causes of death")
    rstpm2::markov_msm(list(model1,smoothpwc), ..., trans = transmat)
}
## use popmort
colon2 = inner_join(survival::colon |> filter(etype==2),
                    mutate(popmort,sex=2-sex,rate) |> filter(year==2000), by=c("age","sex")) |>
    mutate(t=time/365.25)
excess = gsm(Surv(t,status)~factor(rx)+bhazard(rate), data=colon2, df=3)
smoothpwc1 = with(filter(popmort,sex==1 & year==2000),smoothpwc(age+0.5-70,rate)) # example is for men aged 70 years
rs = relative_survival(excess, smoothpwc1, newdata=data.frame(rx="Obs"), t = seq(0,7, length=301))
plot(rs,ggplot=TRUE)




## Spline interpolation
library(schumaker)
d <- data.frame(x=(0:9)+0.5,y=((1:10)-5.5)^2)
xs=seq(0,20,length=301)
plot(xs,Schumaker(d$x,d$y,Extrapolation="Constant")$Spline(xs),type="l")
points(y~x,data=d)
library(splines)
plot(xs,splinefun(d$x,d$y,method="natural")(xs),type="l")
points(y~x,data=d)
d2 <- rbind(d,transform(tail(d,1)[rep(1,2),],x=c(9.5+1e-7,9.5+2e-7)))
xs=seq(9.5-1e-3,9.5+1e-3,length=301)
plot(xs,splinefun(d2$x,d2$y,method="natural")(xs),type="l")
points(y~x,data=d2)
abline(h=20.25,lty=2)
splinefun(d2$x,d2$y,method="natural")(15)-20.25
splinefunx <- function(x, y, method="natural", constant.left=FALSE, constant.right=FALSE, ...) {
    xstar <- x
    ystar <- y
    if (constant.right) {
        xstar <- c(xstar,tail(x,1)+(1:10)*1e-7)
        ystar <- c(ystar,rep(tail(y,1),10))
    }
    if (constant.left) {
        xstar <- c(xstar,x[1]-(1:10)*1e-7)
        ystar <- c(ystar,rep(y[1],10))
    }
    splinefun(xstar, ystar, ..., method=method)
}
splinefunx(d$x,d$y,constant.right=TRUE)(15)-20.25
splinefunx(d$x,d$y,constant.right=TRUE)(15000)-20.25
splinefunx(d$x,d$y)(9.5)-20.25
##
library(numDeriv)
xs=seq(0,20,length=301)
plot(xs, numDeriv::grad(splinefunx(d$x,d$y,constant.right=TRUE), xs), type="l")
rug(d$x)
xs=seq(9.5-1e-3,9.5+1e-3,length=301)
plot(xs, numDeriv::grad(splinefunx(d$x,d$y,constant.right=TRUE), xs), type="l")

## Rather than using mid-points, we could use intervals and then assume a monotone smoother...
numdiff <- function(f,x,eps=1e-5) (f(x+eps)-f(x-eps))/2/eps
d <- transform(data.frame(x0=0:9,x1=(0:9)+1,delta=1),
               xmid=(x0+x1)/2,
               yint=((x1-5)^3-(x0-5)^3)/3)
set.seed(12345)
d <- transform(d,
               yint_ = yint*rgamma(NROW(d),10,10))
xs=seq(0,20,length=301)
## cumd <- transform(d, H=cumsum(c(0,diff(x0)*head(y0,-1))))
cumd <- with(d, data.frame(x=c(x0,tail(x1,1)), H=cumsum(c(0,diff(c(x0,tail(x1,1)))*yint))))
sf.natural <- splinefun(cumd$x,cumd$H, method="natural")
sf.mono <- splinefun(cumd$x,cumd$H, method="hyman")
plot(xs, numdiff(sf.natural,xs), type="l", ylim=c(0,25)) # constant from rhs of last interval
lines(xs, splinefunx(d$xmid,d$yint,constant.right=TRUE)(xs), col="blue") # constant from mid-point
curve((x-5)^2, col="green", add=TRUE)
lines(xs, numdiff(sf.mono,xs), type="l", col="orange") # extrapolates nicely if quadratic

## now with variation in the outcome
cumd <- with(d, data.frame(x=c(x0,tail(x1,1)+1), H=cumsum(c(0,diff(c(x0,tail(x1,1)+1))*yint_))))
sf <- splinefun(cumd$x,cumd$H, method="hyman")
plot(xs, numdiff(sf,xs), type="l", ylim=c(0,35))
lines(xs, splinefunx(d$xmid,d$yint_,constant.right=TRUE)(xs), col="blue")
curve((x-5)^2, col="green", add=TRUE)

## gradient for penalty in the AFT model
library(Rcpp)
library(RcppArmadillo)
fdiff = \(f,x,eps=1e-5)
    sapply(1:length(x),
           \(i) (f("[<-"(x,i,x[i]+eps)) - f("[<-"(x,i,x[i]-eps)))/2/eps)
src = "#include <RcppArmadillo.h>
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat differenceMatrix(int n) {
  return join_rows(zeros(n-1,1),eye(n-1,n-1)) - join_rows(eye(n-1,n-1),zeros(n-1,1));
}
// [[Rcpp::export]]
vec quadraticPenalty(mat Q, vec beta) {
  vec delta = differenceMatrix(Q.n_rows) * Q * beta;
  return delta % delta % (delta<0);
}
// [[Rcpp::export]]
vec gradientPenalty(mat Q, vec beta) { // Q: (nbeta+2) x nbeta 
  size_t n = Q.n_rows;
  mat D = join_rows(zeros(n-1,1),eye(n-1,n-1)) - join_rows(eye(n-1,n-1),zeros(n-1,1)); // (nbeta+1) x (nbeta+2)
  vec delta = D * Q * beta; // nbeta+1
  mat M = Q.t() * D.row(0).t() * D.row(0) * Q * (delta(0)<0.0); // nbeta x nbeta
  for(size_t j=1; j<delta.size(); j++) {
    if (delta(j)<0.0)
      M += Q.t() * D.row(j).t() * D.row(j) * Q;
  }
  return 2*M*beta;
}
// [[Rcpp::export]]
vec gradientPenalty2(mat Q, vec beta) { // Q: (nbeta+2) x nbeta; beta: nbeta
  vec betaStar = Q * beta; // nbeta+2
  size_t n = Q.n_rows; // nbeta+2
  mat D(n-1,n); // nbeta+1 x nbeta+2
  D.zeros();
  for(size_t i=0; i<n-1; i++)
    if (betaStar(i)>betaStar(i+1)) {
      D(i,i) = -1.0;
      D(i,i+1) = 1.0;
    }
  mat M = Q.t() * D.t() * D * Q; // nbeta x nbeta
  return 2*M*beta;
}
"
sourceCpp(code=src)
##
library(rstpm2)
Qmat = attr(nsx(1:10,df=3,intercept = TRUE),"q.const") # 5x3
differenceMatrix(nrow(Qmat)) # 4x5
quadraticPenalty(Qmat, c(1,-1,1)) # 4x1
as.vector(gradientPenalty(Qmat, c(1,-1,1))) 
as.vector(gradientPenalty2(Qmat, c(1,-1,1))) 
fdiff(\(x) sum(quadraticPenalty(Qmat, x)), c(1,-1,1), 1e-5)
as.vector(gradientPenalty(Qmat, c(1,-1,1))) - 
    fdiff(\(x) sum(quadraticPenalty(Qmat, x)), c(1,-1,1), 1e-3)

quadraticPenalty(Qmat, c(1,-1,2)) # 4x1
gradientPenalty(Qmat, c(1,-1,2))

## Error report by Joshua (now fixed)
## Load package
library(rstpm2)
## Define start time
brcancer2 <- transform(brcancer,
                       startTime = ifelse(hormon == 0, rectime * 0.5, 0))
## Stpm2 call without explicit specification of the type argument
stpm2(Surv(startTime, rectime, censrec == 1) ~ hormon,
      data = brcancer2,
      df = 3)
## Stpm2 call with explicit specification of the type argument
## Do not run, this will cause your R-session to terminate
## debug(gsm)
stpm2(Surv(startTime, rectime, censrec == 1, type = "counting") ~ hormon,
      data = brcancer2,
      df = 3)

## Are 2D thin-plate spline approximations in mgcv dense matrices?
library(mgcv)
df=expand.grid(u=seq(0,1,length=100),v=seq(0,1,length=100))
set.seed(12345)
df=transform(df,y=rnorm(nrow(df),sin(u*2*pi)*cos(v*2*pi),0.1))
gam1=gam(y~s(u,v,bs="ts"),data=df)

## Email from Grace
library(rstpm2)
fit = gsm(Surv(rectime,censrec==1)~hormon,data=brcancer,df=3)
fit@lm
attr(fit@lm$terms, "predvars")
nd = data.frame(rectime=1000, hormon=0)
rstpm2:::lpmatrix.lm(fit@lm, nd) %*% coef(fit) # or, more simply
predict(fit@lm, nd)
##
eps=1e-4
(predict(fit@lm, data.frame(rectime=1000+eps, hormon=0)) -
 predict(fit@lm, data.frame(rectime=1000-eps, hormon=0)))/2/eps
## or 
(predict(fit@lm, data.frame(rectime=1000*exp(eps), hormon=0)) -
 predict(fit@lm, data.frame(rectime=1000*exp(-eps), hormon=0)))/2/eps/1000
##
exp(predict(fit@lm, nd))*
    (predict(fit@lm, data.frame(rectime=1000*exp(eps), hormon=0)) -
     predict(fit@lm, data.frame(rectime=1000*exp(-eps), hormon=0)))/2/eps/1000
predict(fit, newdata=nd, type="haz")
plot(fit,type="haz", newdata=data.frame(hormon=0))

## biostat3: coxphHaz and coxphHazList 
as.data.frame.coxphHaz = function(x, row.names=NULL, optional = FALSE, ...) {
    newdata = attr(x,"newdata")
    ## To avoid "row names were found from a short variable and have been discarded":
    rownames(newdata) = NULL 
    data.frame(newdata, x=x$x, y=x$y)
}
as.data.frame.coxphHazList = function(x, row.names=NULL, optional = FALSE, ...) {
    do.call(rbind, lapply(x, as.data.frame))
}
fit <- coxph(Surv(surv_mm/12,status=="Dead: cancer")~agegrp, data=colon)
newdata <- data.frame(agegrp=levels(colon$agegrp))
haz <- suppressWarnings(coxphHaz(fit,newdata))
library(ggplot2)
head(as.data.frame(haz))
ggplot(as.data.frame(haz), aes(x=x,y=y,ymin=y.lower,ymax=y.upper,col=agegrp,fill=agegrp)) +
    geom_ribbon(alpha=0.5) + geom_line(col="black") # + facet_grid(~agegrp)



## Constrained MLE



## relative survival and LEL
library(rstpm2)
library(biostat3)
head(melanoma)
head(popmort)
merged = merge(transform(melanoma, ye = floor(yexit), se=ifelse(sex=="Male",1,2), event=status != "Alive"),
               popmort, by.x=c("age","ye","se"), by.y=c("_age","_year","sex"), all.x=TRUE) 
merged = subset(merged, sex="Female")
## sum(is.na(merged$rate))
model1 = stpm2(Surv(surv_mm,event)~bhazard(rate),data=merged,df=3) # only intercept - was buggy
## model2 = stpm2(Surv(surv_mm,event)~age+bhazard(rate)+cluster(age),data=merged,df=3)
## model1 = stpm2(Surv(surv_mm,event)~1,data=merged,df=3,bhazard=merged$rate)

## mean(subset(merged,event)$yexit)
## We need event counts and person-time in the popmort data-frame
library(dplyr)
temp = cbind(filter(popmort,sex==1) %>% mutate(males=rate,rate=NULL,sex=NULL),
             data.frame(females=filter(popmort,sex==2)$rate))
plot(males~females,temp,log="xy")
abline(0,1,col=2)

##
d = data.frame(y=c(1000,2000),sex=0:1,pt=1e5) # some data
fit = glm(y~sex+offset(log(pt)),data=d, family=poisson) # Poisson regression
pred = predict(fit,newdata=transform(d,pt=1),type="link",se.fit=TRUE)
with(pred,
     exp(data.frame(Estimate=fit,Lower=fit-1.96*se.fit,Upper=fit+1-96*se.fit)))
library(rstpm2)
exp(confint(predictnl(fit,predict,newdata=transform(d,pt=1),type="response")))

set.seed(12345)
n = 1000
x = factor(rbinom(n,1,0.5))
y = rnorm(n, x==1)
fit = lm(y~x)
predict(fit, newdata=data.frame(x=0)) # fails
predict(fit, newdata=data.frame(x=1)) # fails
predict(fit, newdata=data.frame(x=factor(0:1))) # okay
predict(fit, newdata=data.frame(x=factor(0,0:1))) # okay

##
library(expm)
##
K = 2
Z = c(1,0,0)
Y = c(1000,2000); N = 1e4; Sigma=diag(Y/N^2)
lambda=Y/N
##
alpha=cbind(c(-sum(lambda),lambda),matrix(0,K+1,K))
expm(alpha*10) # ok
##
Z = c(1,rep(0,(1+2)*3-1))
alpha=cbind(c(-sum(lambda),lambda,-1,-1,1,0,0,1),matrix(0,length(Z),length(Z)-1))
expm(alpha)
m=matrix(expm(alpha)[4:9,1],2,3)
t(m) %*% Sigma %*% m # includes covariance terms
expm(alpha*10)
m=matrix(expm(alpha*10)[4:9,1],2,3)
sqrt(diag(t(m) %*% Sigma %*% m)) # includes covariance terms


## This approach works!
## We can add further differential equations for the covariance matrix for P (for a time-dependent covariance matrix for the parameters)
library(deSolve)
func <- function(t, y, parms, ...) {
    ## browser()
    index = findInterval(t,parms$t)
    lambda = parms$lambda[index,] # length=K
    K = length(lambda)
    q = rbind(c(-sum(lambda),lambda),matrix(0,K,K+1)) # length=K+1
    Qm = lapply(1:K,function(k) {
        m = matrix(0,K+1,K+1)
        m[1,1] = -1
        m[1,k+1] = 1
        m
    })
    Sigma = parms$Sigma[[index]] # dim=K^2
    P = y[1:(K+1)] # K+1 rowvec
    Pm = matrix(y[(2+K):(1+K+K*(K+1))],K+1,K) # (K+1)*K
    dPdt = P %*% q
    dPmdt = sapply(1:K, function(k) Pm[,k] %*% q + P %*% Qm[[k]]) # (K+1)*K
    dphidt = dPmdt %*% (Sigma %*% t(Pm)) + Pm %*% (Sigma %*% t(dPmdt)) # (K+1)^2
    list(c(dPdt, dPmdt, dphidt)) # length=1+K+(K+1)*K+(K+1)^2
}
Y = c(1000,2000); N = 1e4; Sigma=diag(Y/N^2); K=2
lambda = Y/N
y = c(1,rep(0,K+K*(K+1)+(K+1)^2))
parms = list(t=c(0,Inf), lambda=rbind(lambda), Sigma=list(Sigma))
ode1 = ode(y,seq(0,10),func,parms)
ode1[11,2:4]
sqrt(diag(matrix(ode1[11,11:19],3,3)))
ode1

## Using markov_msm (with an *identity* link for comparability)
fit1=glm(rate~1, data=data.frame(rate=0.1,N=1e4), family=poisson(link="identity"), weight=N)
fit2=glm(rate~1, data=data.frame(rate=0.2,N=1e4), family=poisson(link="identity"), weight=N)
test=markov_msm(list(fit1,fit2),
                trans=matrix(c(NA,1,2,NA,NA,NA,NA,NA,NA),3,3,TRUE), t=0:10, newdata=data.frame(N=1),
                tmvar="t")
print(test,se=TRUE)
test$Pu

## fit1=glm(rate~factor(findInterval(t,0:2))-1, data=data.frame(t=0:2,rate=c(0.1,0.2,0.3),N=1e4), family=poisson(link="identity"), weight=N)

fit1=glm(Y~offset(log(N)), data=data.frame(Y=1000,N=1e4), family=poisson)
fit2=glm(Y~offset(log(N)), data=data.frame(Y=2000,N=1e4), family=poisson)
test=markov_msm(list(fit1,fit2),
                trans=matrix(c(NA,1,2,NA,NA,NA,NA,NA,NA),3,3,TRUE), t=0:10, newdata=data.frame(N=1),
                tmvar="t")
print(test,se=TRUE)


## Inquiry
library(survival)
library(rstpm2)
m <- rstpm2::stpm2(Surv(time, status) ~ sex, data = lung)
plot(m, newdata=data.frame(sex=2))
predict(m, newdata = data.frame(sex = 2, time = 100), type = "rmst", se.fit=TRUE)
## plot(m, newdata=data.frame(sex=2), type="rmst") # fails
## predict(m, newdata = data.frame(sex = 1:2, time = 100), type = "rmst", se.fit=TRUE) # fails

## Inquiry
## Cross-validation with bias correction (Davison and Hinkley, 1997, p. 295)
crossValidation = function(object,data,prediction,K = 10,correction=TRUE) {
    N = nrow(data)
    ndex = 1:N
    part = as.factor(sample(1:K, N, replace = TRUE))
    folds = split(ndex, part)
    p = sapply(folds, length)/N
    sum(sapply(1:K, function(k) {
        training = data[unlist(folds[-k]), ]
        validation = data[folds[[k]], ]
        objectStar = update(object,data=training)
        prediction(objectStar,validation) - p[k]*prediction(objectStar,data)*correction
    })) + prediction(object,data)*correction
}
library(rstpm2)
setMethod("update", "stpm2", function(object, ...) { # I have added this to rstpm2 on GitHub
    object@call = object@Call
    update.default(object, ...)
})
fit = gsm(Surv(rectime,censrec==1)~hormon+x1+x2+x3,data=brcancer,df=3)
set.seed(420)
## Currently, there is no simpler way to predict out-of-sample -2*log-likelihood for an stpm2 model
predictDeviance = function(eventVar) function(object,newdata) {
    events = newdata[[eventVar]]>0
    -2*(sum(predict(object,newdata=newdata[events,,drop=FALSE],
                    type="loghazard"))-
        sum(predict(object,newdata=newdata,type="cumhaz")))
}
crossValidation(fit, brcancer,predictDeviance("censrec"))
set.seed(420)
crossValidation(fit, brcancer,predictDeviance("censrec"),correction=FALSE)
AIC(fit)
BIC(fit)

library(rstpm2)
library(survival)
library(timereg)
library(ggplot2)
library(lattice)
## Two states: Initial -> Final
## Note: this shows how to use markov_msm to estimate survival and risk probabilities based on
## smooth hazard models.
two_states <- function(model, ...) {
    transmat = matrix(c(NA,1,NA,NA),2,2,byrow=TRUE)
    rownames(transmat) <- colnames(transmat) <- c("Initial","Final")
    markov_sde(list(model), ..., trans = transmat)
}
## Note: the first argument is the hazard model. The other arguments are arguments to the
## markov_msm function, except for the transition matrix, which is defined by the new function.
colon2 <- transform(survival::colon, Obs=(rx=="Obs")+0, Lev=(rx=="Lev")+0,Lev_5FU=(rx=="Lev+5FU")+0)
death = aalen(Surv(time,status)~Lev+Lev_5FU, data=subset(colon2,etype==2))
cr = two_states(death, newdata=data.frame(Lev=0,Lev_5FU=0))
plot(cr,ggplot=TRUE,stacked=FALSE,which="L") # ok
plot(cr,ggplot=TRUE,which="L") # ok - is this sensible?
plot(cr) # ok -- no legend
plot(cr,ggplot=TRUE) # ok
plot(cr,lattice=TRUE) # ok -- no legend
plot(cr,ggplot=TRUE,stacked=FALSE) # ok
##
competing_risks <- function(models, ...) {
    transmat = matrix(c(NA,1,2,
                        NA,NA,NA,
                        NA,NA,NA),3,3,byrow=TRUE)
    rownames(transmat) <- colnames(transmat) <- c("Initial","Event 1","Event 2")
    markov_sde(models, ..., trans = transmat)
}
## Note: the first argument is the hazard model. The other arguments are arguments to the
## markov_msm function, except for the transition matrix, which is defined by the new function.
colon2 <- transform(survival::colon, Obs=(rx=="Obs")+0, Lev=(rx=="Lev")+0,Lev_5FU=(rx=="Lev+5FU")+0)
progression = aalen(Surv(time,status)~Lev+Lev_5FU, data=subset(colon2,etype==1))
death = aalen(Surv(time,status)~Lev+Lev_5FU, data=subset(colon2,etype==2))
cr = competing_risks(list(progression,death), newdata=data.frame(Lev=0,Lev_5FU=0))
plot(cr)
##
illness_death <- function(models, ...) {
    transmat = matrix(c(NA,1,2,
                        NA,NA,3,
                        NA,NA,NA),3,3,byrow=TRUE)
    rownames(transmat) <- colnames(transmat) <- c("Initial","Illness","Death")
    markov_sde(models, ..., trans = transmat)
}
## Note: the first argument is the hazard model. The other arguments are arguments to the
## markov_msm function, except for the transition matrix, which is defined by the new function.
colon2 <- transform(survival::colon, Obs=(rx=="Obs"), Lev=(rx=="Lev"),Lev_5FU=(rx=="Lev+5FU"))
d1 = subset(colon2,etype==1)
d2 = subset(colon2,etype==2)
index = d1$time < d2$time & d1$status==1 # index for initial -> recurrence
d12 = d1 # initial -> recurrence
d23 = transform(d2[index,], enter=d1$time[index]) # recurrence -> death
d13 = transform(d2, # initial -> death
                time=pmin(d1$time,d2$time),
                status=ifelse(d1$time==d2$time,d2$status,0))
progression = aalen(Surv(time,status)~Lev+Lev_5FU, data=d12)
death1 = aalen(Surv(time,status)~Lev+Lev_5FU, data=d13)
death2 = aalen(Surv(enter,time,status)~Lev+Lev_5FU, data=d23)
cr = illness_death(list(progression,death1,death2), newdata=data.frame(Lev=FALSE,Lev_5FU=FALSE)) # fails
plot(cr)
##
illness_death = function(models, ...) {
    transmat = matrix(c(NA,1,2,NA,
                        NA,NA,NA,3,
                        NA,NA,NA,NA,
                        NA,NA,NA,NA),4,4,byrow=TRUE)
    rownames(transmat) <- colnames(transmat) <- c("Initial","Illness","DirectDeath","DeathAfterRecurrence")
    markov_sde(models, ..., trans = transmat)
}
## Note: the first argument is the hazard model. The other arguments are arguments to the
## markov_msm function, except for the transition matrix, which is defined by the new function.
colon2 = transform(survival::colon, Obs=(rx=="Obs")+0, Lev=(rx=="Lev")+0, Lev_5FU=(rx=="Lev+5FU")+0)
d1 = subset(colon2,etype==1)
d2 = subset(colon2,etype==2)
index = d1$time < d2$time & d1$status==1 # index for initial -> recurrence
d12 = d1 # initial -> recurrence
d23 = transform(d2[index,], enter=d1$time[index]) # recurrence -> death
d13 = transform(d2, # initial -> death
                time=pmin(d1$time,d2$time),
                status=ifelse(d1$time==d2$time,d2$status,0))
progression = aalen(Surv(time,status)~Lev+Lev_5FU, data=d12)
death1 = aalen(Surv(time,status)~Lev+Lev_5FU, data=d13)
death2 = aalen(Surv(enter,time,status)~Lev+Lev_5FU, data=d23)
cr = illness_death(list(progression,death1,death2), newdata=data.frame(Lev=0,Lev_5FU=0:1),
                   weights=c(-1,1), los=FALSE)
plot(cr,ggplot=TRUE) + facet_grid(~Lev_5FU) # ok
plot(cr,ggplot=TRUE,stacked=FALSE) + facet_grid(state~Lev_5FU) # ok
plot(cr,ggplot=TRUE,stacked=FALSE,which="L") + facet_grid(state~Lev_5FU) # ok
plot(standardise(cr),stacked=FALSE,ggplot=TRUE) # CIs are *very* wide -- error?

## Non-parametric baseline: SDE approach due to Ryalen and colleagues
markov_sde <- function(models, trans, newdata, init=NULL, nLebesgue=1e4+1, los=FALSE, nOut=300,
                        weights=1) {
    transfun <- function(tmat) {
        indices <- sort(as.vector(tmat)); indices <- setdiff(indices,NA)
        nStates <- nrow(tmat)
        out <- do.call(rbind,
                lapply(indices, function(i) {
                    index2 <- which(tmat == i)
                    from <- (index2-1) %% nStates +1
                    to <- (index2-1) %/% nStates + 1
                    data.frame(from=from,to=to)
                }))
        matrix(as.integer(as.matrix(out)),nrow(out))-1L
    }
    ## TODO check parameters
    nStates <- nrow(trans)
    nTrans <- sum(!is.na(trans))
    if (is.null(init)) {
        init <- c(1,rep(0,nStates-1))
    }
    stopifnot(length(init)==nStates)
    ## init <- rep(init,nrow(newdata))
    n <- sum(sapply(models,attr,"orig.max.clust"))
    cumHazList <- lapply(models, function(object)
        -t(log(predict(object, newdata=newdata, se=FALSE)$S0)))
    timesList <- lapply(models, function(model) model$cum[,1])
    eventTimes <- times <- sort(unique(unlist(timesList)))
    hazList <- lapply(1:nrow(newdata), function(i) {
        hazMatrix <- matrix(0,nrow=nTrans,ncol=length(times))
        for (j in 1:nTrans) {
            hazMatrix[j,match(timesList[[j]],times)] <- diff(c(0,cumHazList[[j]][,i]))
        }
        hazMatrix
    })
    hazMatrix <- do.call(rbind,hazList)
    if (length(weights)==1 && weights != 0)
        weights <- rep(weights,nrow(newdata))/(weights*nrow(newdata))
    vcov <- matrix(1,1,1)
    if (!los) {
        out <- .Call("plugin_P_by", n, nrow(newdata), hazMatrix, init, transfun(trans), weights,
                     vcov,
                     PACKAGE="rstpm2")
        out$P <- out$X
        out$P.se <- sqrt(out$variance)
    } else {
        out <- .Call("plugin_P_L_by",
                     n, nrow(newdata), hazMatrix, init, transfun(trans), times, weights, nOut,
                     vcov, nLebesgue,
                     PACKAGE="rstpm2")
        PIndex <- 1:(nrow(out$X)/2)
        tr <- function(x) array(as.vector(x), dim=c(nStates,nrow(newdata),ncol(x)))
        out$P <- tr(out$X[PIndex,])
        out$L <- tr(out$X[-PIndex,])
        out$P.se <- sqrt(tr(out$variance[PIndex,]))
        out$L.se <- sqrt(tr(out$variance[-PIndex,]))
    }
    out$n <- n
    out$times <- if (los) out$time else times
    out$newdata <- newdata
    out$trans <- trans
    out$los <- los
    out$init <- init
    out$weights <- weights
    class(out) <- "markov_sde"
    if (!all(weights==0)) {
        stand <- out # Warning: copy! This may be a bad idea...
        if (!los) {
            stand$P <- stand$Y
            stand$P.se <- sqrt(stand$varY)
        } else {
            PIndex <- 1:(nrow(out$Y)/2)
            stand$P <- stand$Y[PIndex,]
            stand$L <- stand$Y[-PIndex,]
            stand$P.se <- sqrt(stand$varY[PIndex,])
            stand$L.se <- sqrt(stand$varY[-PIndex,])
        }
        ## tidy up
        stand$X <- stand$Y
        stand$variance <- stand$varY
        ## out$X <- out$variance <- out$Y <- out$varY <- stand$Y <- stand$varY <- NULL
        stand$newdata <- stand$newdata[1,,drop=FALSE]
        class(stand) <- "markov_sde"
        out$stand <- stand
    }
    out
}
standardise.markov_sde <- function(object) {
    structure(object$stand, class="markov_sde")
}
cr = two_states(death, newdata=data.frame(rx=levels(survival::colon$rx)))

cr = two_states(death, newdata=data.frame(Obs=1))


plot.markov_sde <- function(x, y, stacked=TRUE, which=c("P","L"), 
                            xlab="Time", ylab=NULL, col=2:6, border=col,
                            ggplot2=FALSE, lattice=FALSE, alpha=0.2,
                            strata=NULL,
                            ...) {
    stopifnot(inherits(x,"markov_sde"))
    which <- match.arg(which)
    if (!missing(y)) warning("y argument is ignored")
    ## ylab defaults
    if (is.null(ylab))
        ylab <- if(which=='P') "Probability" else "Length of stay"
    if (ggplot2)
        rstpm2:::ggplot.markov_msm(x, which=which, stacked=stacked, xlab=xlab, ylab=ylab,
                          alpha=alpha, ...)
    else if (lattice)
        rstpm2:::xyplot.markov_msm(x, which=which, stacked=stacked, xlab=xlab, ylab=ylab,
                          col=col, border=border, strata=strata, ...)
    else {
        if (is.null(index) && nrow(x$newdata)>1) {
            warning("More than one set of covariates; defaults to weighted estimator")
            x <- x$stand # Warning: replacement
            index <- 1
        }
        browser()
        if (is.null(index)) index <- 1
        df <- merge(x$newdata[index,,drop=FALSE], as.data.frame(x))
        states <- unique(df$state)
        if (stacked) {
            out <- graphics::plot(range(x$times, na.rm=TRUE),0:1, type="n", xlab=xlab, ylab=ylab, ...)
            lower <- 0
            for (i in length(states):1) { # put the last state at the bottom
                df2 <- df[df$state==states[i],]
                if (length(lower)==1) lower <- rep(0,nrow(df2))
                upper <- lower+df2[[which]]
                graphics::polygon(c(df2$time,rev(df2$time)), c(lower,rev(upper)),
                                  border=border[i], col=col[i])
                lower <- upper
            }
            graphics::box()
            invisible(out)
        }
        else stop('Unstacked plot not implemented in base graphics; use ggplot2=TRUE or lattice=TRUE')
    }
}

plot(standardise(cr)) # ok (no legend)
plot(standardise(cr), ggplot=TRUE) # ok (includes legend)
plot(standardise(cr), lattice=TRUE) # ok (no legend)
plot(cr, ggplot=TRUE) + facet_grid(~ rx) # ok

plot(cr, index=1) # incorrect
plot(cr, ggplot=TRUE, index=1) # incorrect
plot(cr, lattice=TRUE, index=1) # incorrect
plot(standardise(cr), ggplot=TRUE, stacked=FALSE) # incorrect: wrong CIs for standardised estimates
plot(cr, ggplot=TRUE, stacked=FALSE) + facet_grid(~ rx) # incorrect


as.data.frame.markov_sde <- function(x, row.names=NULL, ci=TRUE,
                                      P.conf.type="logit", L.conf.type="log",
                                      P.range=c(0,1), L.range=c(0,Inf),
                                      ...) {
    if (any(x$weights<0))
        P.conf.type <- L.conf.type <- "plain"
    .id. <- 1:nrow(x$newdata)
    nStates <- nrow(x$trans)
    state.names <- rownames(x$trans)
    stateNames <- if (!is.null(rownames(x$trans))) rownames(x$trans) else 1:nrow(x$trans)
    out <- expand.grid(state=stateNames, .id.=.id., time=x$times)
    out <- cbind(x$newdata[out$.id.,],out)
    names(out)[1:ncol(x$newdata)] <- colnames(x$newdata)
    out$P <- as.vector(x$P)
    out$P.se <- as.vector(x$P.se)
    out <- out[order(out$.id.,out$state,out$time),]
    out$.id. <- NULL
    if (ci) {
        tmp <- rstpm2:::surv.confint(out$P,out$P.se, conf.type=P.conf.type, min.value=P.range[1], max.value=P.range[2])
        out$P.lower <- tmp$lower
        out$P.upper <- tmp$upper
    }
    if (x$los) {
        out$L <- as.vector(x$L)
        out$L.se <- as.vector(x$L.se)
        if (ci) {
            tmp <- rstpm2:::surv.confint(out$L,out$L.se, conf.type=L.conf.type, min.value=L.range[1], max.value=L.range[2])
            out$L.lower <- tmp$lower
            out$L.upper <- tmp$upper
        }
    }
    if(!is.null(row.names)) rownames(out) <- row.names
    out
}
temp = as.data.frame(cr)
head(temp)
ggplot(temp,aes(x=time,y=P,fill=rx,ymin=P.lower,ymax=P.upper)) + geom_line() +
    geom_ribbon(alpha=0.5) + facet_grid(state~rx) # ok
rstpm2:::ggplot.markov_msm(temp) + facet_grid(~rx) # ok


## Competing risks
## Note: this shows how to adapt the markov_msm model for competing risks.
competing_risks <- function(listOfModels, ...) {
    nRisks = length(listOfModels)
    transmat = matrix(NA,nRisks+1,nRisks+1)
    transmat[1,1+(1:nRisks)] = 1:nRisks
    rownames(transmat) <- colnames(transmat) <- c("Initial",names(listOfModels))
    rstpm2::markov_msm(listOfModels, ..., trans = transmat)
}
## Note: The first argument for competing_risks is a list of models. Names from that list are
## used for labelling the states. The other arguments are as per the markov_msm function,
## except for the transition matrix, which is defined by the competing_risks function.
recurrence = gsm(Surv(time,status)~factor(rx), data=survival::colon, subset=(etype==1), df=3)
death = gsm(Surv(time,status)~factor(rx), data=survival::colon, subset=(etype==2), df=3)
cr = competing_risks(list(Recurrence=recurrence,Death=death),
                     newdata=data.frame(rx=levels(survival::colon$rx)),
                     t = seq(0,2500, length=301))
## Plot the probabilities for each state for three different treatment arms
plot(cr, ggplot=TRUE) + facet_grid(~ rx)
## And: differences in probabilities
cr_diff = diff(subset(cr,rx=="Lev+5FU"),subset(cr,rx=="Obs"))
plot(cr_diff, ggplot=TRUE, stacked=FALSE)


## Example using Crowther and Lambert (2018)
library(readstata13)
library(transform.hazards)
library(timereg)
library(rstpm2) # standardise
mex.1 <- read.dta13("https://fmwww.bc.edu/repec/bocode/m/multistate_example.dta")
transmat <- rbind("Post-surgery"=c(NA,1,2), 
                  "Relapsed"=c(NA,NA,3),
                  "Died"=c(NA,NA,NA))
colnames(transmat) <- rownames(transmat)
mex.2 <- transform(mex.1,osi=(osi=="deceased")+0)
levels(mex.2$size)[2] <- ">20-50 mm" # fix typo
mex <- mstate::msprep(time=c(NA,"rf","os"),status=c(NA,"rfi","osi"),
                      data=mex.2,trans=transmat,id="pid",
                      keep=c("age","size","nodes","pr_1","hormon"))
mex <- transform(mex,
                 size2=(unclass(size)==2)+0, # avoids issues with TRUE/FALSE
                 size3=(unclass(size)==3)+0,
                 hormon=(hormon=="yes")+0,
                 Tstart=Tstart/12,
                 Tstop=Tstop/12)
## Slow fitting...
c.ar <- aalen(Surv(Tstart,Tstop,status) ~ const(age) + size2 + size3 + const(nodes) + pr_1 + const(hormon),
              data = subset(mex, trans==1))
c.ad <- aalen(Surv(Tstart, Tstop, status) ~ const(age) + const(size) + const(nodes) + const(pr_1) + const(hormon),
              data = subset(mex, trans==2))
c.rd <- aalen( Surv(Tstart,Tstop,status) ~ const(age) + const(size) + const(nodes) + pr_1 + const(hormon),
              data=subset(mex, trans==3))
##
nd <- expand.grid(nodes=seq(0,20,10), size=levels(mex$size))
nd <- transform(nd, age=54, pr_1=3, hormon=0, size2=(unclass(size)==2)+0, size3=(unclass(size)==3)+0)
system.time(fit1 <- markov_sde(list(c.ar,c.ad,c.rd), trans=transmat, newdata=nd[c(1,2),], los=TRUE))
system.time(fit0 <- markov_sde(list(c.ar,c.ad,c.rd), trans=transmat, newdata=nd[c(1,2),]))
plot(fit1)
plot(fit1, ggplot=TRUE)
plot(fit1, lattice=TRUE, stacked=FALSE, which="L")
## plot(fit1,which="L",xlim=NULL) # does this make sense?
df1 <- as.data.frame(fit1) 
fit2 <- standardise(fit1)
plot(fit2,ggplot=TRUE)
df2 <- as.data.frame(fit2)
diff1 <- diff(fit1)
plot(diff1)


## type="af"
library(rstpm2)
fit = gsm(Surv(rectime,censrec==1)~hormon,data=brcancer,df=3)
plot(fit,type="af",newdata=brcancer,
    exposed=function(data) transform(data,hormon=1))
pred2 <- predict(fit,type="af",newdata=brcancer, grid=TRUE,
    exposed=function(data) transform(data,hormon=1),
    se.fit=TRUE,full=TRUE)
with(pred2, matplot(rectime,cbind(Estimate,lower,upper),type="l",lty=c(1,2,2), col=1,
                    xlab="Time since treatment (days)", ylab="PAF", ylim=c(0,0.5)))

## Bug report from Joshua
# Loading rstpm2 package
library(rstpm2)
# Fit stmp2 model using breastcancer data set
fpm_model <- stpm2(Surv(rectime, censrec) ~ hormon,
                   data = brcancer,
                   df = 3,
                   tvc = list("hormon" = 3))
# Predict hazard difference comparing hormon users and non-users
# using full=TRUE option for obtaining a full data set for ggplot()
predict(
  fpm_model,
  type = "hdiff",
  newdata = data.frame(hormon = 0),
  var = "hormon",
  grid = TRUE,
  se.fit = TRUE,
  full = TRUE)

## test plots for PO models with random effects
library(rstpm2)
set.seed(12345)
logit <- binomial()$linkfun
expit <- binomial()$linkinv
Spo <- function(eta) expit(eta)
eta <- function(t) {}
rPOgamma <- function(n,eta,theta,eps=1e-16) {
    ## U <- pmin(runif(n),1-eps)
    U <- runif(n)
    ## solve_t(U=expit(-eta(t))^theta), where eta(0)=-Inf and eta(Inf)=Inf
    V <- -logit(U^(1/theta))
    if (any(is.na(V))) browser()
    c(list(U=U),vuniroot(function(t) eta(t)-V, lower=rep(1e-50,n), upper=rep(1e100,n)))
}
t0 <- rPOgamma(1e5, function(t) log(t), 1)
with(t0,plot(U,root,log="y"))
t1 <- rPOgamma(n=1e5, function(t) log(t), rgamma(10,1))
range(t0$root)
range(t1$root)
plot(density(log(t0$root)))

## Plots with three levels
library(rstpm2)
table(brcancer$x4) # cancer stage
## define indicators for the "exposed" levels (*R* needs these to be numeric)
d <- transform(brcancer, x4.2=(x4==2)+0, x4.3=(x4==3)+0)
## fit the model with the indicators as main effects and with tvc
fit <- stpm2(Surv(rectime,censrec==1)~x4.2+x4.3,data=d,df=3,tvc=list(x4.2=2,x4.3=2))
## predict for each exposure level
pred2 <- predict(fit,newdata=data.frame(x4.2=0,x4.3=0),type="hr",var="x4.2", grid=TRUE, full=TRUE,
                 se.fit=TRUE)
pred3 <- predict(fit,newdata=data.frame(x4.2=0,x4.3=0),type="hr",var="x4.3", grid=TRUE, full=TRUE,
                 se.fit=TRUE)
pred <- transform(rbind(pred2,pred3), x4=factor(ifelse(x4.2,2,3)))
library(ggplot2)
ggplot(pred, aes(x=rectime, y=Estimate, ymin=lower, ymax=upper, fill=x4, col=x4)) +
    geom_line() + coord_cartesian(ylim=c(0,10)) + geom_ribbon(alpha=0.3, col=NA)



## how to extract the design information
library(rstpm2)
fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,df=3)
names(fit@model.frame)
attributes(fit@model.frame[[3]])
fit@lm$terms
fit@x <- matrix()
ls(environment(fit@model.frame))
##
withEnvs <-  sapply(slotNames(fit),function(nm) !is.null(environment(slot(fit,nm))))
slotNames(fit)[withEnvs]
lapply(slotNames(fit)[withEnvs],function(nm) ls(environment(slot(fit,nm))))
##
withEnvs <-  sapply(fit@args,function(obj) !is.null(environment(obj)))
names(fit@args)[withEnvs]
lapply(fit@args[withEnvs],function(obj) ls(environment(obj)))


## bug in predict for meansurv
library(rstpm2)
library(ggplot2)
fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,df=3)
## Easy to use plot and lines functions
plot(fit, newdata=transform(brcancer, hormon=0), type="meansurv")
lines(fit, newdata=transform(brcancer, hormon=1), type="meansurv",lty=2)
## More tedious to do for different covariate patterns with ggplot2
pred0 <- predict(fit, newdata=transform(brcancer, hormon=0), type="meansurv", full=TRUE, se.fit=TRUE,
                 grid=TRUE)
pred1 <- predict(fit, newdata=transform(brcancer, hormon=1), type="meansurv", full=TRUE, se.fit=TRUE,
                 grid=TRUE)
## bug with values returned AsIs - I'll try to fix this
unAsIs <- function(object)
    if (inherits(object,"AsIs")) "class<-"(object, setdiff(class(object), "AsIs")) else object
pred <- rbind(transform(unAsIs(pred1),hormon=1),transform(unAsIs(pred0),hormon=0))
pred <- transform(pred, Hormone=ifelse(hormon==1,"Yes","No"))
ggplot(pred, aes(x=rectime,y=Estimate,ymin=lower,ymax=upper,fill=Hormone)) +
    xlab("Time since diagnosis (years)") +
    ylab("Standardised survival") +
    geom_ribbon(alpha=0.6) +
    geom_line()


## bug in predict for meansurv
library(devtools)
install_github("mclements/rstpm2", ref="develop")
library(rstpm2)
fit.tvc <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,df=3, tvc=list(hormon=3))
out <- predict(fit.tvc, newdata=transform(brcancer,hormon=1),type="meansurv",grid=TRUE, se.fit=TRUE,
               full=TRUE)

## bug in plot(..., xlab="Something")
library(rstpm2)
fit <- stpm2(Surv(rectime, censrec)~hormon, data=brcancer,df=3)
plot(fit, newdata=data.frame(hormon=1), type="hr", xlab="Time since diagnosis (years)", var="hormon",
     ylab="Hazard ratio", main="Lung cancer-BMI")

## predict linear predictor
library(rstpm2)
fit <- stpm2(Surv(rectime, censrec)~hormon, data=brcancer,df=3)
predict(fit, newdata=data.frame(hormon=0:1, rectime=1000), type="link")
predict(fit, newdata=data.frame(hormon=1, rectime=1000), type="link") -
    predict(fit, newdata=data.frame(hormon=0, rectime=1000), type="link")

##
library(rstpm2)
library(Hmisc)
d <- local({
    set.seed(12345)
    x <- rep(0:1,length=200)
    y <- rexp(length(x), exp(-3+x))
    data.frame(x,y,e=TRUE)
    })
fit0 <- stpm2(Surv(y, e)~1, data=d,df=3)
fit <- stpm2(Surv(y, e)~x, data=d,df=3)
rcorr.cens(predict(fit0,type="link")-predict(fit0,newdata=transform(d,x=0),type="link"),
           with(d,Surv(y,e)))
rcorr.cens(-predict(fit,type="link")+
           predict(fit,newdata=transform(d,x=0),type="link"),
           with(d,Surv(y,e)))
summary(fit)

## testing - Bug requires loading devtools *before* rstpm2
setwd("~/src/R/rstpm2")
library(devtools)
devtools::test()

## offset
library(rstpm2)
brcancer2 <- transform(brcancer,off=0.1)
fit <- stpm2(Surv(rectime,censrec==1)~hormon, df = 2, data=brcancer2)
fit2 <- stpm2(Surv(rectime,censrec==1)~hormon+offset(off), df = 2, data=brcancer2)
fit
fit2
zeroModel(fit)

## p-value for survival differences
library(rstpm2)
fit <- stpm2(Surv(rectime,censrec==1)~hormon, df = 2, data=brcancer,tvc=list(hormon=2))
test <- predict(fit, type="sdiff", var="hormon", newdata=data.frame(hormon=0,rectime=500),se.fit=TRUE)
z <- test[1,1]/((test[1,3]-test[1,2])/2/qnorm(0.975))
2*pnorm(-abs(z))
##
## test the difference in survival rates
test <- predictnl(fit, function(object,newdata=NULL) {
  lp1 <- predict(object, newdata=data.frame(hormon=1,rectime=500), type="surv")
  lp2 <- predict(object, newdata=data.frame(hormon=0,rectime=500), type="surv")
  lp1-lp2
  })
with(test, c(fit=fit,
             se.fit=se.fit,
             statistic=fit/se.fit,
             p=2*pnorm(-abs(fit/se.fit))))
## same p-value as predict(..., type="sdiff")
##
s18 <- summary(survfit(Surv(rectime,censrec==1)~hormon,data=brcancer),time=500)
z2<-diff(s18$surv)/sqrt(sum(s18$std.err^2))
2*pnorm(-abs(z2))
##
library(bpcp)
with(brcancer,
     fixtdiff(rectime,censrec,hormon,500,doall = TRUE))

## bug #3 on GitHub
library(rstpm2)
# Function with model_formula paramter **ERRORS**
f_bad <- function(model_formula)
  stpm2(formula = model_formula, df = 2, data=brcancer, link.type = "PH")
# Function with formula paramter to match the name pf stpm2's formula parameter **PASSES**
f_good <- function(formula)
  stpm2(formula = formula, df = 2, data=brcancer, link.type = "PH")
f_bad (model_formula = Surv(rectime,censrec==1)~hormon)
f_good (formula = Surv(rectime,censrec==1)~hormon)


## random draws
library(rstpm2)
predict.cumhaz <-
          function(object, newdata=NULL, ...)
{
    stopifnot(inherits(object,"stpm2") || inherits(object,"pstpm2"))
    args <- object@args
    beta <- coef(object)
    if (is.null(newdata))
        X <- args$X
    else if (inherits(object, "stpm2")) {
          X <- object@args$transX(lpmatrix.lm(object@lm, newdata), newdata)
      }
    else if (inherits(object, "pstpm2")) {
           X <- object@args$transX(predict(object@gam, newdata, type="lpmatrix"), newdata)
      }
    link <- object@link # cf. link for transformation of the predictions
    eta <- as.vector(X %*% beta)
    link$H(eta)
}
simulate <- function(object, nsim=nrow(as.data.frame(newdata)),
                     newdata=as.data.frame(object@data), lower=1e-6, upper=1e5, ...) {
    e <- rexp(nsim)
    objective <- function(time) {
        newdata[[object@timeVar]] <- time
        predict.cumhaz(object, type="cumhaz", newdata=newdata) - e
    }
    vuniroot(objective, lower=rep(lower,length=nsim), upper=rep(upper,length=nsim))$root
}
fit1 <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer)

set.seed(12345)
d <- do.call(rbind,lapply(1:100,function(i) brcancer))
system.time(r <- simulate(fit1, newdata=d))
length(r)

library(rstpm2)
fit1 <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer)
negll2 <- function(beta,svalues,missingIndicator,k=2)
    sapply(svalues, function(s)
        fit1@args$logli2(coef(fit1),ifelse(missingIndicator, s*coef(fit1)[k], 0)))
head(negll2(coef(fit1),c(0,1),(1:686) %% 2))


## values for the tests
library(rstpm2)
brcancer2 <- transform(rstpm2::brcancer, w=1)
brcancer2$w[1] <- NA
fit1 <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer2,weights=w)
nd <- data.frame(hormon=1,rectime=1000)
predict(fit1, newdata=nd, type="surv")
predict(fit1, newdata=nd, type="fail")
predict(fit1, newdata=nd, type="haz")
predict(fit1, newdata=nd, type="hr",var="hormon")
predict(fit1, newdata=nd, type="hdiff", var="hormon")


## predictions for relative survival (email from Anke Richters)
set.seed(12345)
d <- with(list(t0=rexp(10000), # constant hazard
               t1=rweibull(10000, 1.5), # rising hazard
               bg=rexp(2*10000)), # constant background hazard (rate=1)
          data.frame(t=pmin(bg,c(t0,t1)), x=rep(0:1,each=10000), e=TRUE, bhazard=1))
library(rstpm2)
uniroot(function(x) pexp(x)-pweibull(x,shape=1.5), c(1e-6,10)) # survival overlaps at 1
uniroot(function(x) 1-dweibull(x,shape=1.5)/pweibull(x,shape=1.5,lower=FALSE),
        c(1e-6,10)) # hazards overlap at 0.444
fit <- stpm2(Surv(t,e)~x,data=d,tvc=list(x=3),bhazard=d$bhazard)
## fit <- pstpm2(Surv(t,e)~1,data=d,smooth.formula=~s(log(t))+s(log(t),by=x),bhazard=d$bhazard)
ts <- seq(0,4,length=301)[-1]
## hazard plot - ok
plot(fit,newdata=data.frame(x=0),type="hazard",ylim=c(0,4))
lines(fit,newdata=data.frame(x=1),type="hazard",lty=2)
abline(h=1,col="blue") # theorical
lines(ts,dweibull(ts,shape=1.5)/pweibull(ts,shape=1.5,lower.tail=FALSE),lty=2,col="blue") # theorical
## survival plot - ok
plot(fit,newdata=data.frame(x=0),type="surv")
lines(fit,newdata=data.frame(x=1),type="surv",lty=2)
lines(ts, pexp(ts,lower.tail=FALSE), col="blue") # theoretical
lines(ts, pweibull(ts,shape=1.5,lower.tail=FALSE), lty=2, col="blue") # theoretical
## hr - now fixed
plot(fit,newdata=data.frame(x=0),type="hr",var="x")
lines(ts, dweibull(ts,shape=1.5)/pweibull(ts,shape=1.5,lower.tail=FALSE), lty=2, col="blue") # theoretical
abline(h=1,lty=2)
## plot(fit,newdata=data.frame(x=0),type="hr",exposed=function(data) transform(data,x=1)) # same mistake
plot(fit,newdata=data.frame(x=0),type="sdiff",var="x",ylim=c(-1,1))
S0 <- predict(fit,newdata=data.frame(x=0),grid=TRUE,full=TRUE)
S1 <- predict(fit,newdata=data.frame(x=1),grid=TRUE,full=TRUE)
lines(S0$t,S1$Estimate-S0$Estimate,col="blue")
abline(h=0)
abline(v=1)

library(rstpm2)
fit1 <- coxph(Surv(rectime,censrec==1)~ns(x1,df=2),data=brcancer)


## Missing values in predictions
library(rstpm2)
brcancer2 <- brcancer
brcancer$rectime[1] <- NA
fit1 <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer2)
pred1 <- predict(fit1,newdata=brcancer2)
head(pred1)
fit1
summary(fit1)
brcancer2 <- brcancer
brcancer$hormon[1] <- NA
fit1 <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer2)
pred1 <- predict(fit1,newdata=brcancer2)
head(pred1)
fit1
summary(fit1)
##
library(rstpm2)
brcancer2 <- brcancer
brcancer$rectime[1] <- NA
fit1 <- pstpm2(Surv(rectime,censrec==1)~hormon,data=brcancer2)
pred1 <- predict(fit1,newdata=brcancer2)
head(pred1)
fit1
summary(fit1)
brcancer2 <- brcancer
brcancer$hormon[1] <- NA
fit1 <- pstpm2(Surv(rectime,censrec==1)~hormon,data=brcancer2)
pred1 <- predict(fit1,newdata=brcancer2)
head(pred1)
fit1
summary(fit1)


## Examples using predictnl for Alessandro
library(rstpm2)
brcancer2 <- transform(brcancer, x4.23=x4 %in% 2:3)
fit1 <- stpm2(Surv(rectime,censrec==1)~hormon*x4.23,data=brcancer2,df=3)
summary(fit1)
newd <- data.frame(hormon=0,x4.23=FALSE)
plot(fit1, newdata=newd)
RERI <- function(object, newdata,
                 var1, val1=1, 
                 var2, val2=1) {
    exp1 <- function(data) {data[[var1]] <- val1; data}
    exp2 <- function(data) {data[[var2]] <- val2; data}
    s00 <- predict(object, newdata, type="surv")
    s10 <- predict(object, newdata=exp1(newdata), type="surv")
    s01 <- predict(object, newdata=exp2(newdata), type="surv")
    s11 <- predict(object, newdata=exp1(exp2(newdata)), type="surv")
    -(s11-s10-s01+s00)/(1-s00)
}
times <- seq(0,2500,length=301)[-1]
reri <- RERI(fit1,newdata=transform(newd,rectime=times),var1="hormon",var2="x4.23",val2=TRUE)
plot(times,reri,type="l")
reri2 <- predictnl(fit1,fun=RERI,newdata=transform(newd,rectime=times),var1="hormon",var2="x4.23",val2=TRUE)
with(reri2, matplot(times,fit+cbind(0,-1.96*se.fit,+1.96*se.fit),type="l",lty=c(1,2,2),col=1,
                    xlab="Time since diagnosis", ylab="RERI"))
abline(h=0,lty=3)

RERI.hr <- function(object, newdata,
                 var1, val1=1, 
                 var2, val2=1) {
    exp1 <- function(data) {data[[var1]] <- data[[var1]]+val1; data}
    exp2 <- function(data) {data[[var2]] <- data[[var2]]+val2; data}
    h00 <- predict(object, newdata, type="haz")
    h10 <- predict(object, newdata=exp1(newdata), type="haz")
    h01 <- predict(object, newdata=exp2(newdata), type="haz")
    h11 <- predict(object, newdata=exp1(exp2(newdata)), type="haz")
    (h11-h10-h01+h00)/h00
}
RERI.hr(fit1,newdata=transform(newd,rectime=1000),var1="hormon",var2="x4.23",val2=TRUE)
predictnl(fit1,fun=RERI.hr,newdata=transform(newd,rectime=1000),var1="hormon",var2="x4.23",val2=TRUE)



## testing of relative survival
library(rstpm2)
ayear <- 365.24
brcancer2 <- transform(brcancer, age=80*ayear, sex="male", year=as.Date("1980-01-01"), time=1, recyear=rectime/ayear)
rate0 <- survexp(time~1,data=brcancer2,method="individual.h",scale=ayear)
(fit1 <- stpm2(Surv(recyear,censrec==1)~hormon,data=brcancer2,df=2,cure=T,bhazard=rate0))
head(predict(fit1,type.relsurv="excess"))
head(predict(fit1,type.relsurv="total"))
head(brcancer2)

ayear <- 365.24
timeVar <- substitute(times)
scale <- ayear
rmap <- substitute(list())
newdata <- data.frame(sex=c("male",rep("male",5)),age=ayear*60,year=2002,times=c(1,1:5))
survexp1 <- do.call(survexp, list(substitute(I(timeVar*scale)~1,list(timeVar=timeVar)),
                                  ratetable=survexp.us,
                                  scale=scale,
                                  rmap=rmap,
                                  cohort=FALSE,
                                  data=newdata))


plot(fit1, newdata=data.frame(hormon=1,age=80,sex="male",year=1980))
## lines(fit1, newdata=data.frame(hormon=1,age=80,sex="male",year=1980))

## Bug report from Alessandro for 1.4.0
library(rstpm2)
data(kidney)
fitg = stpm2(Surv(time, status) ~ age + sex, cluster = kidney$id, data = kidney, 
  RandDist = "Gamma")
head(predict(fitg))
fitln = stpm2(Surv(time, status) ~ age + sex, cluster = kidney$id, data = kidney, 
  RandDist = "LogN")
head(predict(fitln))
fitln = stpm2(Surv(time, status) ~ age + sex, cluster = kidney$id, data = kidney, Z=~age-1,
  RandDist = "LogN")
head(predict(fitln))

## test meanhr
library(rstpm2)
fit <- stpm2(Surv(rectime, censrec==1) ~ x4+x5, data = brcancer, df=3)
fit <- stpm2(Surv(rectime, censrec==1) ~ x4+x5, data = brcancer, df=3)
summary(fit)
eform(fit)
plot(fit, newdata=data.frame(hormon=0,x4=0,x5=0))
plot(fit, newdata=data.frame(hormon=0,x4=0,x5=0),type="hazard")
plot(fit, newdata=data.frame(hormon=0,x4=0,x5=0), type="hr", exposed=function(data) transform(data, x4=1))
plot(fit, newdata=transform(brcancer,x4=1), type="meanhr", exposed=function(data) transform(data, x4=2))
plot(fit, newdata=transform(brcancer,x4=1), type="meanhaz")

## test rmst
library(rstpm2)
fit <- stpm2(Surv(rectime, censrec==1) ~ hormon, data = brcancer, df=3)
plot(fit, newdata=data.frame(hormon=1))
predict(fit, newdata=data.frame(hormon=1,rectime=1000), type="rmst", se.fit=TRUE)
predict(fit, newdata=data.frame(hormon=0,rectime=1000), type="rmst", se.fit=TRUE)

library(devtools)
install.packages("bbmle")
devtools::install_github("mclements/rstpm2",ref="develop")

## 2017-06-21 
## Verify: the choice of basis dimension (default: k=10) for penalized regression splines is not sensitive to estimates
## Adjusted by a constant coefficient (e.g. alpha=2) to correct potential overfitting by GCV for lambda
## alpha = 1.4 suggested by Kim and Gu (2004)
library(rstpm2)
## k = 7
pfit7 <- pstpm2(Surv(rectime, censrec==1) ~ hormon, data = brcancer, smooth.formula = ~ s(log(rectime), k=7), alpha=2)
plot(pfit7, newdata = data.frame(hormon=0), type="hazard")

## k = 27
pfit27 <- pstpm2(Surv(rectime, censrec==1) ~ hormon, data = brcancer, smooth.formula = ~ s(log(rectime), k=27), alpha=2)
plot(pfit27, newdata = data.frame(hormon=0), type="hazard")

## Estimated effective degree of freedom (EDF)
pfit7@edf  ## 5.36
pfit27@edf ## 5.96


require(coxme)

## Fix error in code for gradli
library(rstpm2)
data(brcancer)
fit <- stpm2(Surv(rectime,censrec) ~ hormon,data=transform(brcancer,censrec=1))
fit <- stpm2(Surv(rectime,censrec==1) ~ hormon,data=brcancer,cure=TRUE)
fit <- stpm2(Surv(rectime,censrec==1) ~ hormon,data=brcancer)

plot(fit,newdata=data.frame(hormon=1),type="uncured",exposed=function(data) transform(data,rectime=2500))
X <- fit@args$X
XD <- fit@args$XD
args <- fit@args
beta.est <- coef(fit)
eta <- as.vector(X %*% beta.est)
etaD <- as.vector(XD %*% beta.est)
link <- switch(fit@args$link,PH=rstpm2:::link.PH,PO=rstpm2:::link.PO)
h <- link$h(eta,etaD)  # - as.vector(predict(fit, type="haz")) ## Ok!
H <- link$H(eta)  #- as.vector(predict(fit, type="cumhaz")) ## Ok!
gradh <- as.matrix(link$gradh(eta,etaD, args))
gradH <- as.matrix(link$gradH(eta, args))

gradli <- residuals(fit, type="gradli") ## n*npar
dim(gradli)
gradli2 <- gradH - ifelse(fit@args$event,1/h,0)*gradh
head(gradli + gradli2)


## Gamma frailty
require(rstpm2)
brcancer2 <- transform(brcancer, id=rep(1:(nrow(brcancer)/2),each=2))
brcancer2$hormon[1] <- NA
fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer2, cluster=brcancer2$id)
summary(fit)
plot(fit,newdata=data.frame(hormon=1),type="margsurv")
fit2 <- pstpm2(Surv(rectime,censrec==1)~hormon,data=brcancer2, cluster=brcancer2$id)
summary(fit2)
plot(fit2,newdata=data.frame(hormon=1),type="margsurv")

# Aranda-Ordaz link
refresh
require(rstpm2)
## PH
summary(fit <- stpm2(Surv(rectime,censrec==1)~1,data=brcancer,link="PH", df=3))
summary(fit <- stpm2(Surv(rectime,censrec==1)~1,data=brcancer,link="AO", df=3)) # Same: OK
summary(fit <- stpm2(Surv(rectime,censrec==1)~1,data=brcancer,link="PO", df=3))
summary(fit <- stpm2(Surv(rectime,censrec==1)~1,data=brcancer,link="AO", theta.AO=1, df=3)) # Same: OK
summary(fit <- pstpm2(Surv(rectime,censrec==1)~1,data=brcancer,link="AO", theta.AO=0.5))

refresh
require(rstpm2)
## PH
summary(fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,link="PH", df=3))
predict(fit, newdata=transform(brcancer,rectime=1000),type="meansurv",keep.attributes=FALSE,se.fit=TRUE,use.gr=F)
predict(fit, newdata=transform(brcancer,rectime=1000),type="meansurv",keep.attributes=FALSE,se.fit=TRUE,use.gr=T)
plot(fit,newdata=transform(brcancer,hormon=1),type="meansurv",times=seq(10,1500,by=10))
plot(fit,newdata=transform(brcancer,hormon=2),type="meansurv",times=seq(10,1500,by=10),lty=2,add=TRUE)

newd <- merge(transform(brcancer,rectime=NULL), data.frame(rectime=c(500,1000)))
unlist(predict(fit,newdata=newd,type="af",exposed=function(data) transform(data,hormon=1),keep.attributes=FALSE,se.fit=TRUE) -
       predict(fit,newdata=newd,type="af",exposed=function(data) transform(data,hormon=1),keep.attributes=FALSE,se.fit=TRUE,use.gr=FALSE))

system.time(plot(fit,type="af",exposed=function(data) transform(data,hormon=1),recent=TRUE))
system.time(plot(fit,type="af",exposed=function(data) transform(data,hormon=1),recent=FALSE))
plot(fit,newdata=NULL,type="meansurv",ci=F)
plot(fit,newdata=NULL,type="meansurvdiff",exposed=function(data) transform(data,hormon=1))

plot(fit,newdata=data.frame(hormon=1),type="surv",ci=F)
plot(fit,newdata=data.frame(hormon=1),type="fail",ci=T)

unlist(predict(fit,newdata=newd,type="meansurvdiff",exposed=function(data) transform(data,hormon=1),keep.attributes=FALSE,se.fit=TRUE) -
      predict(fit,newdata=newd,type="meansurvdiff",exposed=function(data) transform(data,hormon=1),keep.attributes=FALSE,se.fit=TRUE,use.gr=FALSE))
unlist(predict(fit,newdata=newd,type="meansurv",keep.attributes=FALSE,se.fit=TRUE)-
       predict(fit,newdata=newd,type="meansurv",keep.attributes=FALSE,se.fit=TRUE,use.gr=FALSE))


## comparison with AF
## Example 1: clustered data with frailty U
require(AF)
set.seed(12345)
expit <- function(x) 1 / (1 + exp( - x))
n <- 100
m <- 2
alpha <- 1.5
eta <- 1
phi <- 0.5
beta <- 1
id <- rep(1:n,each=m)
U <- rep(rgamma(n, shape = 1 / phi, scale = phi), each = m)
Z <- rnorm(n * m)
X <- rbinom(n * m, size = 1, prob = expit(Z))
## Reparametrize scale as in rweibull function
weibull.scale <- alpha / (U * exp(beta * X)) ^ (1 / eta)
t <- rweibull(n * m, shape = eta, scale = weibull.scale)
## Right censoring
cen <- runif(n * m, 0, 10)
delta <- as.numeric(t < cen)
t <- pmin(t, cen)
d <- data.frame(t, delta, X, Z, id)


require(rstpm2)
fit2 <- stpm2(formula = Surv(t, delta) ~ X + Z + X * Z, data = d, df=1, cluster=d$id, smooth.formula=~log(t))
predict(fit2, type="af", newdata=transform(d,t=1),exposed=function(data) transform(data, X=0), se.fit=TRUE)
plot(fit2, type="af", exposed=function(data) transform(data, X=0))
plot(fit2, type="meansurvdiff", exposed=function(data) transform(data, X=0))
plot(fit2, type="meansurv")

fit3 <- pstpm2(formula = Surv(t, delta) ~ X + Z + X * Z, data = d, df=1, cluster=d$id)
predict(fit3, type="af", newdata=transform(d,t=1),exposed=function(data) transform(data, X=0), se.fit=TRUE)
plot(fit3, type="meansurv")

predict(fit2, newdata=transform(d,t=1), type="meansurv")

## check analytical gradients for margsurv and marghaz
predict(fit2, newdata=data.frame(t=1,X=1,Z=1), type="margsurv", use.gr=TRUE, se.fit=TRUE)-predict(fit2, newdata=data.frame(t=1,X=1,Z=1), type="margsurv", use.gr=FALSE, se.fit=TRUE)
predict(fit2, newdata=data.frame(t=1,X=1,Z=1), type="marghaz", use.gr=TRUE, se.fit=TRUE)-predict(fit2, newdata=data.frame(t=1,X=1,Z=1), type="marghaz", use.gr=FALSE, se.fit=TRUE)
predict(fit2, newdata=data.frame(t=1,X=1,Z=1), type="hazard", use.gr=TRUE, se.fit=TRUE)-predict(fit2, newdata=data.frame(t=1,X=1,Z=1), type="hazard", use.gr=FALSE, se.fit=TRUE)

require(boot)
meansurv <- function(data,index) predict(fit2, newdata=transform(data[index,,drop=FALSE],t=1), type="meansurv")
meansurv(d,TRUE)
boot1 <- boot(d, meansurv, R=1000)
boot.ci(boot1)

require(rstpm2)
fit <- stpm2(formula = Surv(t, delta) ~ X + Z + X * Z, data = d, df=1)
diag(vcov(fit))
fit <- stpm2(formula = Surv(t, delta) ~ X + Z + X * Z, data = d, frailty=FALSE, cluster=d$id, df=1)
diag(vcov(fit))
fit <- stpm2(formula = Surv(t, delta) ~ X + Z + X * Z, data = d, cluster=d$id, df=1)
diag(vcov(fit))
##
fit <- stpm2(formula = Surv(t, delta) ~ X + Z + X * Z, data = d, cluster = d$id, df=1)
predict(fit,type="af",newdata=transform(d,t=1),exposed=function(data) transform(data,X=0),keep.attributes=FALSE,se.fit=TRUE) 
fit <- stpm2(formula = Surv(t, delta) ~ X + Z + X * Z, data = d, cluster = d$id, df=4)
predict(fit,type="af",newdata=transform(d,t=1),exposed=function(data) transform(data,X=0),keep.attributes=FALSE,se.fit=TRUE)

fit <- stpm2(formula = Surv(t, delta) ~ X + Z + X * Z, data = d, cluster=d$id, df=1)
predict(fit,type="af",newdata=transform(d,t=1),exposed=function(data) transform(data,X=0),keep.attributes=FALSE,se.fit=TRUE) 

## Fit a frailty object
library(stdReg)
fit <- stdReg::parfrailty(formula = Surv(t, delta) ~ X + Z + X * Z, data = d, clusterid = "id")
summary(fit)
## Estimate the attributable fraction from the fitted frailty model
time <- c(seq(from = 0.2, to = 1, by = 0.2))
time <- 1
## debug(AFfrailty)
library(AF)
AFfrailty_est <- AFparfrailty(object = fit, data = d, exposure = "X", times = time, clusterid = "id")
AFfrailty_est
##AF:::summary.AF(AFfrailty_est)





## tvc for Maarten Coemans
require(rstpm2)
brcancer <- transform(brcancer, x1c=x1-mean(x1))
summary(fit.tvc <- stpm2(Surv(rectime,censrec==1)~ hormon+x1c, data=brcancer, df=3,  tvc=list(hormon=2,x1c=2)))
plot(fit.tvc,newdata=data.frame(hormon=0,x1c=-10),type="hr", var="hormon")
plot(fit.tvc,newdata=data.frame(hormon=0,x1c=+10),type="hr", var="hormon", add=TRUE,ci=FALSE,line.col=2)
## same model
summary(fit.tvc <- stpm2(Surv(rectime,censrec==1)~ hormon+x1c, data=brcancer,
                         smooth.formula = ~ns(log(rectime),df=3)+hormon:ns(log(rectime),df=2)+x1c:ns(log(rectime),df=2)))
## and again...
summary(fit.tvc <- stpm2(Surv(rectime,censrec==1)~ x1c, data=brcancer,
                         smooth.formula = ~ns(log(rectime),df=3)+hormon:ns(log(rectime),df=3,intercept=TRUE)+x1c:ns(log(rectime),df=2)))
## new model with time different time transformation for the TVCs
summary(fit.tvc <- stpm2(Surv(rectime,censrec==1)~ hormon+x1c, data=brcancer,
                         smooth.formula = ~ns(log(rectime),df=3)+hormon:ns(rectime,df=2)+x1c:ns(rectime,df=2)))
plot(fit.tvc,newdata=data.frame(hormon=0,x1c=-10),type="hr", var="hormon")
plot(fit.tvc,newdata=data.frame(hormon=0,x1c=+10),type="hr", var="hormon", add=TRUE,ci=FALSE,line.col=2)
##
## not including the main effect and no intercept is not the same
summary(fit.tvc <- stpm2(Surv(rectime,censrec==1)~ x1c, data=brcancer,
                         smooth.formula = ~ns(log(rectime),df=3)+hormon:ns(log(rectime),df=2)+x1c:ns(log(rectime),df=2)))

## Standardised survival
require(rstpm2)
plot.meansurv <- function(x, y=NULL, times=NULL, newdata=NULL, add=FALSE, ci=!add, rug=!add, recent=FALSE,
                          xlab=NULL, ylab="Mean survival", lty=1, line.col=1, ci.col="grey", ...) {
    if (is.null(times)) stop("plot.meansurv: times argument should be specified")
    if (is.null(newdata)) newdata <- x@data
    times <- times[times !=0]
    if (recent) {
        newdata <- do.call("rbind",
                           lapply(times, 
                                  function(time) {
                                      newdata[[x@timeVar]] <- newdata[[x@timeVar]]*0+time
                                      newdata
                                  }))
        pred <- predict(x, newdata=newdata, type="meansurv", se.fit=ci) # requires recent version
        pred <- if (ci) rbind(c(Estimate=1,lower=1,upper=1),pred) else c(1,pred)
    } else {
        pred <- lapply(times, 
                       function(time) {
                           newdata[[x@timeVar]] <- newdata[[x@timeVar]]*0+time
                           predict(x, newdata=newdata, type="meansurv", se.fit=ci)
                       })
        pred <- if (ci) rbind(c(Estimate=1,lower=1,upper=1),do.call("rbind", pred)) else c(1,unlist(pred))
        }
    times <- c(0,times)
    if (is.null(xlab)) xlab <- deparse(x@timeExpr)
    if (!add) matplot(times, pred, type="n", xlab=xlab, ylab=ylab, ...)
    if (ci) {
        polygon(c(times,rev(times)),c(pred$lower,rev(pred$upper)),col=ci.col,border=ci.col)
        lines(times,pred$Estimate,col=line.col,lty=lty)
    } else {
        lines(times,pred,col=line.col,lty=lty)
    }
    if (rug) {
        Y <- x@y
        eventTimes <- Y[Y[,ncol(Y)]==1,ncol(Y)-1]
        rug(eventTimes,col=line.col)
    }
    return(invisible(y))
}
brcancer <- transform(brcancer, x1c=x1-mean(x1))
summary(fit.tvc <- stpm2(Surv(rectime,censrec==1)~ hormon+x1c, data=brcancer, df=3,  tvc=list(hormon=2,x1c=2)))
times <- seq(0,3000,by=100)
plot.meansurv(fit.tvc, newdata=transform(brcancer, hormon=1), times=times,
              ylim=c(0.2,1))
plot.meansurv(fit.tvc, times=times, newdata=transform(brcancer, hormon=0), line.col=2, add=TRUE)


## Examples using ns() for covariates - this was buggy.
refresh
require(rstpm2)
summary(fit <- stpm2(Surv(rectime,censrec==1)~1,
                     smooth.formula=~ns(log(rectime),df=3)+ns(x1,df=3),
                     data=brcancer,link="PH"))
summary(fit <- stpm2(Surv(rectime,censrec==1)~ns(x1,df=3), df=3,data=brcancer,link="PH"))
summary(fit <- pstpm2(Surv(rectime,censrec==1)~ns(x1,df=3), data=brcancer,link="PH"))


grad <- function(f,x,eps=1e-5)
    sapply(1:length(x), function(i) {
        lower <- upper <- x
        upper[i] <- x[i]+eps
        lower[i] <- x[i]-eps
        (f(upper)-f(lower))/2/eps
    })
link <- function(S,theta=0.5) log((S^(-theta)-1)/theta)
S <- ilink <- function(eta,theta=0.5) exp(-log(theta*exp(eta)+1)/theta)
H <- function(eta,theta=0.5) -log(S(eta,theta))
h <- function(eta,etaD,theta=0.5) exp(eta)*etaD/(theta*exp(eta)+1)
gradH <- function(eta,X,theta=0.5) exp(eta)*X/(1+theta*exp(eta))
gradh <- function(eta,etaD,X,XD,theta=0.5) {
    eta <- as.vector(eta)
    etaD <- as.vector(etaD)
    ((theta*exp(2*eta)+exp(eta))*XD+exp(eta)*etaD*X) /
        (theta*exp(eta)+1)^2
}


X <- cbind(1,1:2,1) # (constant, t, x)
XD <- cbind(0,1:2,0)
beta <- c(0.1, 0.2, 0.3)
eta <- as.vector(X %*% beta)
etaD <- as.vector(XD %*% beta)
S(eta)
H(eta)
h(eta,etaD) - grad(function(t) H(cbind(1,t,1) %*% beta), 1) # OK
gradH(eta,X) - grad(function(beta) H(X %*% beta), beta) # OK
gradh(eta,etaD,X,XD)
grad(function(beta) h(X %*% beta, XD %*% beta), beta)

ilink(link(.1))
link(ilink(.1))




require(abind)
X <- matrix(seq(0,1,length=5*10),nrow=10)
beta <- seq(0,1,length=5)
H <- exp(as.vector(X %*% beta))
dHdbeta <- X * H # row=indiv, col=beta

d2Hdbeta2 <- aperm(abind(lapply(1:ncol(X), function(k) X[,k] * X * H),along=3),c(2,3,1))
abind(lapply(1:nrow(X), function(i) (X[i,] %*% t(X[i,])) * H[i]),along=3) -
    aperm(abind(lapply(1:ncol(X), function(k) X[,k] * X * H),along=3),c(2,3,1))

numder <- function(f,x,eps=1e-8) (f(x+eps)-f(x-eps))/2/eps
expit <- function(x) 1/(1+exp(-x))
numder(expit,2)
expit(2)*expit(-2)
numder(dnorm,2)
-dnorm(2)*2

require(mgcv)
d <- data.frame(x = seq(0,1,length=100), x2=rnorm(100), y = rnorm(100))
fit <- gam(y~s(x)+s(x2,by=x), data=d)
X <- predict(fit,d,type="lpmatrix")
X0 <- predict(fit,transform(d,x=0),type="lpmatrix")
Xstar <- X-X0
index0 <- rstpm2:::which.dim(Xstar)
lapply(fit$smooth, function(s) {
    which((1:ncol(X) %in% index0)[s$first.para:s$last.para]) # index for S'
})
lapply(fit$smooth, function(s) {
    range(which((1:ncol(X) %in% s$first.para:s$last.para)[index0]))
})
lapply(fit$smooth,"[[","S")
## outline: given a full index=1:n, a reduced index set index0 and a smoother with first.para, last.para and a square matrix S, return a revised first.para', last.para' and matrix S'
## For S': 



refresh
require(rstpm2)
## additive
fit <- stpm2(Surv(rectime,censrec==1)~1,data=brcancer,link="AH",
             smooth.formula=~ns(rectime,df=4)+hormon:ns(rectime,df=3), optimiser="NelderMead")
summary(fit)
fit2 <- stpm2(Surv(rectime,censrec==1)~1,data=brcancer,link="AH",
             smooth.formula=~ns(rectime,df=4)+hormon:ns(rectime,df=3))
summary(fit2)
plot(fit2,newdata=data.frame(hormon=0),type="haz")
plot(fit2,newdata=data.frame(hormon=1),add=TRUE,lty=2,type="haz")

fit <- pstpm2(Surv(rectime,censrec==1)~1,data=brcancer,link="AH",
             smooth.formula=~s(rectime)+s(rectime,by=hormon))
plot(fit,newdata=data.frame(hormon=0),type="haz")
plot(fit,newdata=data.frame(hormon=1),add=TRUE,lty=2,type="haz")

## test robust estimators from penalized models
## most spline coefficients become statistically significant
require(rstpm2)
summary(pstpm2(Surv(rectime/365,censrec==1)~hormon,data=brcancer,robust=FALSE))
summary(pstpm2(Surv(rectime/365,censrec==1)~hormon,data=brcancer,robust=TRUE))

## robust standard errors for clustered data
refresh
require(rstpm2)
brcancer2 <- transform(brcancer,
                       id=rep(1:(nrow(brcancer)/2),each=2))
fit <- stpm2(Surv(rectime,censrec==1)~1,data=brcancer)
summary(fit)
fit <- stpm2(Surv(rectime,censrec==1)~1,data=brcancer, cluster=brcancer2$id, robust=TRUE)
summary(fit)
##
require(rstpm2)
brcancer2 <- transform(brcancer, id=rep(1:(nrow(brcancer)/2),each=2))
fit <- stpm2(Surv(rectime,censrec==1)~1,data=brcancer, cluster=brcancer2$id)
summary(fit)
predict(fit,type="gradli")



## Stata estimated coef for hormon
## PH:     -.3614357
## PO:     -.474102
## Probit: -.2823338
system.time(print( stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer, stata=TRUE)))
system.time(print(pfit <- pstpm2(Surv(rectime,censrec==1)~hormon,smooth.formula=~s(log(rectime))+s(x1),data=brcancer)))
##
system.time(print( stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,type="PO")))
system.time(print(pstpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,type="PO")))
##
system.time(print( stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,type="probit")))
system.time(print(pstpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,type="probit"))) # slow

summary(fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,smooth.formula=~nsx(log(rectime), df=4, stata.stpm2.compatible = TRUE)))

if (FALSE) {
    debug(pstpm2)
    pfit <- pstpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,sp=1)
    ## towards the end of the pstpm2 function...
    sum(diag(solve(optimHess(coef(mle2),negllsp,sp=1)) %*% optimHess(coef(mle2),negll0sp,sp=1)))
    sum(diag(solve(optimHess(coef(mle2),negllsp,sp=fit$sp)) %*% optimHess(coef(mle2),negll0sp,sp=fit$sp)))
    negllsp(coef(mle2),sp=1)
    negll0sp(coef(mle2),sp=1)
}
update.list <- function(list,...) {
    args <- list(...)
    for (name in names(args))
        list[[name]] <- args[[name]]
    list
}

## right censored
## Stata estimated coef for hormon (PH): -.3614357
refresh
require(rstpm2)
summary(fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,
      smooth.formula=~nsx(log(rectime),df=3,stata=TRUE),trace=0))
summary(fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,
      smooth.formula=~nsx(log(rectime),df=3,stata=TRUE),trace=0,optimiser="NelderMead"))
summary(fit2 <- pstpm2(Surv(rectime,censrec==1)~hormon,data=brcancer))
summary(fit2 <- pstpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,optimiser="NelderMead"))

## delayed entry
## Stata estimated coef for hormon (PH): -1.162504
library(rstpm2)
brcancer2 <- transform(brcancer,startTime=ifelse(hormon==0,rectime/2,0))
## brcancer2 <- transform(brcancer,startTime=0.1)
##debug(rstpm2:::meat.stpm2)
summary(fit <- stpm2(Surv(startTime,rectime,censrec==1)~hormon,data=brcancer2,
      smooth.formula=~nsx(log(rectime),df=3,stata=TRUE)))
summary(fit2 <- pstpm2(Surv(startTime,rectime,censrec==1)~hormon,data=brcancer2,control=list(optimiser="NelderMead"))) # OK
summary(fit2 <- pstpm2(Surv(startTime,rectime,censrec==1)~hormon,data=brcancer2))
plot(fit,newdata=data.frame(hormon=1))
lines(fit2,newdata=data.frame(hormon=1),lty=2)
head(predict(fit)) # OK
head(predict(fit,se.fit=TRUE))
## delayed entry and tvc (problems?)
summary(fit <- stpm2(Surv(startTime,rectime,censrec==1)~hormon,data=brcancer2,
                     smooth.formula=~nsx(rectime,df=3)+hormon:nsx(rectime,df=3,stata=TRUE)))
head(predict(fit,se.fit=TRUE)) 
pstpm2(Surv(startTime,rectime,censrec==1)~hormon,data=brcancer2)

## left truncated with clusters
library(rstpm2)
brcancer2 <- transform(brcancer,
                       startTime=ifelse(hormon==0,rectime/2,0),
                       id=rep(1:(nrow(brcancer)/2),each=2))
##debug(stpm2)
summary(fit <- stpm2(Surv(startTime,rectime,censrec==1)~hormon,data=brcancer2,
                     cluster=brcancer2$id, control=list(optimiser="NelderMead"),recurrent=TRUE,
                     smooth.formula=~nsx(log(rectime),df=3,stata=TRUE)))
summary(fit0 <- stpm2(Surv(startTime,rectime,censrec==1)~hormon,data=brcancer2,
                     control=list(optimiser="NelderMead"),
                     smooth.formula=~nsx(log(rectime),df=3,stata=TRUE)))

require(foreign)
require(rstpm2)
stmixed <- read.dta("http://fmwww.bc.edu/repec/bocode/s/stmixed_example2.dta")
stmixed2 <- transform(stmixed, start = ifelse(treat,stime/2,0))
summary(r2 <- stpm2(Surv(stime,event)~treat,data=stmixed2,cluster=stmixed$trial,RandDist="LogN",df=3,Z=~treat-1,adaptive=TRUE,optimiser="NelderMead"))
summary(r2 <- stpm2(Surv(stime,event)~treat,data=stmixed2,cluster=stmixed$trial,RandDist="LogN",df=3,Z=~treat-1,adaptive=TRUE))
##
library(foreign)
library(rstpm2)
stmixed <- read.dta("http://fmwww.bc.edu/repec/bocode/s/stmixed_example2.dta")
summary(r <- stpm2(Surv(stime,event)~treat,data=stmixed,cluster=stmixed$trial,RandDist="LogN",df=3,Z=~treat-1))

## non-adaptive
system.time(print(summary(r2 <- stpm2(Surv(stime,event)~treat,data=stmixed2,cluster=stmixed$trial,RandDist="LogN",df=3,Z=~treat,adaptive=FALSE,nodes=20,optimiser="NelderMead")))) # slow and gradients not close to zero
system.time(print(summary(r2 <- stpm2(Surv(stime,event)~treat,data=stmixed2,cluster=stmixed$trial,RandDist="LogN",df=3,Z=~treat,nodes=20,adaptive=FALSE)))) # gradients close to zero

## random intercept and random slope with 20 nodes
summary(r2 <- stpm2(Surv(stime,event)~treat,data=stmixed2,cluster=stmixed$trial,RandDist="LogN",df=3,Z=~treat,nodes=20,adaptive=FALSE)) 
summary(r2 <- stpm2(Surv(stime,event)~treat,data=stmixed2,cluster=stmixed$trial,RandDist="LogN",df=3,Z=~treat,adaptive=FALSE,nodes=20)) # gradients close to zero

## Simple examples with no random effects and with a random intercept (check: deviances)
summary(r2 <- stpm2(Surv(stime,event)~treat,data=stmixed,df=3))
summary(r2 <- stpm2(Surv(stime,event)~treat,data=stmixed,cluster=stmixed$trial,RandDist="LogN",df=3,Z=~treat-1))


## check modes and sqrttau
args <- r2@args
args$return_type <- "modes"
.Call("model_output", args, package="rstpm2")
args$return_type <- "variances" # fudge
.Call("model_output", args, package="rstpm2")


## check gradients
args <- r2@args
args$return_type <- "gradient"
.Call("model_output", args, package="rstpm2")
fdgrad <- function(obj,eps=1e-6) {
    args <- obj@args
    args$return_type <- "objective"
    sapply(1:length(args$init), function(i) {
        largs <- args
        largs$init[i] <- args$init[i]+eps
        f1 <- .Call("model_output", largs, package="rstpm2")
        largs$init[i] <- args$init[i]-eps
        f2 <- .Call("model_output", largs, package="rstpm2")
        data.frame(f1,f2,gradient=(f1-f2)/2.0/eps)
    })
}
fdgrad(r2,1e-3)


## random intercept and random slope
summary(r2 <- stpm2(Surv(stime,event)~treat,data=stmixed2,cluster=stmixed$trial,RandDist="LogN",df=3,Z=~treat))

summary(stpm2(Surv(start,stime,event)~treat,data=stmixed2))
summary(r2 <- stpm2(Surv(stime,event)~treat,data=stmixed2,cluster=stmixed$trial,RandDist="LogN"))
##
summary(r2 <- pstpm2(Surv(stime,event)~treat,data=stmixed2,cluster=stmixed$trial,RandDist="LogN"))
##
summary(r2 <- stpm2(Surv(stime,event)~treat+factor(trial),data=stmixed2,cluster=stmixed$trial,RandDist="LogN",Z=~treat-1))
##
summary(r2 <- pstpm2(Surv(stime,event)~treat+factor(trial),data=stmixed2,cluster=stmixed$trial,RandDist="LogN",Z=~treat-1))

## gradients for RE parameters
require(expm)
link.PH <- list(link=function(S) log(-log(S)),
                ilink=function(eta) exp(-exp(eta)),
                h=function(eta,etaD) etaD*exp(eta),
                H=function(eta) exp(eta),
                gradh=function(eta,etaD,obj) obj$XD*exp(eta)+obj$X*etaD*exp(eta),
                gradH=function(eta,obj) obj$X*exp(eta))
link.AH <- list(link=function(S) -log(S),
                ilink=function(eta) exp(-eta),
                h=function(eta,etaD) etaD,
                H=function(eta) eta,
                gradh=function(eta,etaD,obj) obj$XD,
                gradH=function(eta,obj) obj$X)
corrtrans <- function(x) (1-exp(-x)) / (1+exp(-x))
l <- function(gamma=c(log(.3),log(.4),0.5)) {
    ## initially assume an additive model
    H <- function(eta) eta
    h <- function(eta,etaD) etaD
    delta <- 1
    eta <- 1
    etaD <- 2
    Z <- c(1,2)
    d <- c(1,1)
    mode <- c(1,2)
    ## define sqrtmat
    var <- exp(gamma[c(1,2)])
    rho <- corrtrans(gamma[3])
    cov <- rho*sqrt(var[1]*var[2])
    Sigma <- matrix(c(var[1],cov,cov,var[2]),2,2)
    ## SqrtSigma <- chol(Sigma)
    SqrtSigma <- sqrtm(as(Sigma,"symmetricMatrix"))
    b <- mode + SqrtSigma %*% d
    Zb <- sum(Z*b)
    l <- delta*log(h(eta+Zb,etaD))-H(eta+Zb)
    l
}
g <- c(log(.3),log(.4),0.5)
## gradient wrt g using finite differences
sapply(1:3, function(i,eps=1e-6) {
    g[i] <- g[i]+eps
    val <- l(g)
    g[i] <- g[i]-2*eps
    (val - l(g))/2/eps
})

bf <- function(gamma=c(log(.3),log(.4),0.5)) {
    Z <- c(1,2)
    d <- c(1,1)
    mode <- c(1,2)
    ## define sqrtmat
    var <- exp(gamma[c(1,2)])
    rho <- corrtrans(gamma[3])
    cov <- rho*sqrt(var[1]*var[2])
    Sigma <- matrix(c(var[1],cov,cov,var[2]),2,2)
    ## SqrtSigma <- chol(Sigma)
    SqrtSigma <- sqrtm(as(Sigma,"symmetricMatrix"))
    b <- mode + SqrtSigma %*% d
    b
}
## gradient for b wrt g using finite differences
(gradbwrtg <- sapply(1:3, function(i,eps=1e-6) {
    g[i] <- g[i]+eps
    val <- bf(g)
    g[i] <- g[i]-2*eps
    (val - bf(g))/2/eps
}))

var <- exp(g[c(1,2)])
rho <- corrtrans(g[3])
cov <- rho*sqrt(var[1]*var[2])
## wrt g[1]
matrix(c(g[1]*var[1],cov*g[1]/2,cov*g[1]/2,0),2,2)

##
A=1; B=0.5; D=2
M=matrix(c(A,B,B,D),2,2)
tau <- A+D
delta <- A*D-B*B
s <- sqrt(delta)
t <- sqrt(tau+2*s)
sqrtM <- matrix(c(A+s,B,B,D+s),2,2)/t
sqrtM %*% sqrtM

lb <- function(b) {
    ## initially assume an additive model
    H <- function(eta) eta
    h <- function(eta,etaD) etaD
    delta <- 1
    eta <- 1
    etaD <- 2
    Z <- c(1,2)
    d <- c(1,1)
    mode <- c(1,2)
    Zb <- sum(Z*b)
    l <- delta*log(h(eta+Zb,etaD))-H(eta+Zb)
    l
}
## gradient wrt b using finite differences
gradlwrtb <- sapply(1:2, function(i,eps=1e-6) {
    b[i] <- b[i]+eps
    val <- lb(b)
    b[i] <- b[i]-2*eps
    (val - lb(b))/2/eps
})
t(gradbwrtg) %*% gradlwrtb


## gradient for a multivariate normal
require(mvtnorm)
fdgrad <- function(f,x, ..., eps=1.0e-6) 
    sapply(1:length(x), function(i) {
        e <- rep(0,length(x))
        e[i] <- 1
        (f(x+e*eps, ...)-f(x-e*eps, ...))/2/eps
    })
## hess(i,i) = (-f2 +16.0*f1 - 30.0*f0 + 16.0*fm1 -fm2)/(12.0*hi*hi);
fdhessian <- function(f,x, ..., eps=1.0e-5) 
    sapply(1:length(x), function(i) {
        ei <- rep(0,length(x))
        ei[i] <- 1
        sapply(1:length(x), function(j) {
            ej <- rep(0,length(x))
            ej[j] <- 1
            if (i==j) (-f(x+2*ei*eps, ...)+16*f(x+ei*eps, ...)-30*f(x,...)+16*f(x-ei*eps, ...)-f(x-2*ei*eps, ...))/12/eps/eps else (f(x+eps*ei+eps*ej)-f(x+eps*ei-eps*ej)-f(x-eps*ei+eps*ej)+f(x-eps*ei-eps*ej))/4/eps/eps
        })
    })
-fdgrad(mvtnorm::dmvnorm, c(1,1), c(0,0), Sigma <- matrix(c(1,.5,.5,1),2), log=TRUE)
## fdhessian(mvtnorm::dmvnorm, c(1,1), c(0,0), Sigma <- matrix(c(1,.5,.5,1),2))
as.vector(solve(Sigma) %*% c(1,1))


summary(fit <- stpm2(Surv(start,stime,event)~treat,data=stmixed2,cluster=stmixed$trial,optimiser="NelderMead",recurrent=TRUE))
summary(fit <- stpm2(Surv(start,stime,event)~treat,data=stmixed2,cluster=stmixed$trial,recurrent=TRUE))
summary(fit <- stpm2(Surv(start,stime,event)~treat,data=stmixed2,cluster=stmixed$trial,RandDist="LogN",
                   optimiser="NelderMead", recurrent=TRUE))
summary(r2 <- stpm2(Surv(start,stime,event)~treat,data=stmixed2,cluster=stmixed$trial,RandDist="LogN",
                    recurrent=TRUE))

summary(fit <- stpm2(Surv(start,stime,event)~treat,data=stmixed2,cluster=stmixed$trial,optimiser="NelderMead"))
summary(fit <- stpm2(Surv(start,stime,event)~treat,data=stmixed2,cluster=stmixed$trial))
summary(fit <- stpm2(Surv(start,stime,event)~treat,data=stmixed2,cluster=stmixed$trial,RandDist="LogN",
                   optimiser="NelderMead"))
summary(r2 <- stpm2(Surv(start,stime,event)~treat,data=stmixed2,cluster=stmixed$trial,RandDist="LogN"))



summary(r <- stpm2(Surv(start,stime,event)~treat,data=stmixed2))



## check gradients
args <- fit@args
args$return_type <- "gradient"
.Call("model_output", args, package="rstpm2")
fdgrad <- function(obj,eps=1e-6) {
    args <- obj@args
    args$return_type <- "objective"
    sapply(1:length(args$init), function(i) {
        largs <- args
        largs$init[i] <- args$init[i]+eps
        f1 <- .Call("model_output", largs, package="rstpm2")
        largs$init[i] <- args$init[i]-eps
        f2 <- .Call("model_output", largs, package="rstpm2")
        data.frame(f1,f2,gradient=(f1-f2)/2.0/eps)
    })
}
fdgrad(fit,1e-3)



require(rstpm2)
require(mgcv)
x=seq(0,1,length=5001)
set.seed(12345)
y=rnorm(length(x),sin(2*pi*x))
i <- x>0.65
d=data.frame(x=x[i],y=y[i])
fit <- gam(y~s(x),data=d)
## plot(fit)
plot(x,predict(fit,newdata=data.frame(x=x)),type="l")
plot(x,y)

## weighted estimates
refresh
require(rstpm2)
## unequal weights
brcancer2 <- transform(brcancer,w=ifelse(hormon==0,10,1))
## unweighted 
summary(fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer2,
      smooth.formula=~nsx(log(rectime),df=3,stata=TRUE)))
## weighted estimates
## stpm2
summary(stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer2,weights=w,
              smooth.formula=~nsx(log(rectime),df=3,stata=TRUE)))
## stpm2 robust
summary(stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer2,weights=w,robust=TRUE,
      smooth.formula=~nsx(log(rectime),df=3,stata=TRUE)))
summary(pstpm2(Surv(rectime,censrec==1)~hormon,data=brcancer2))
## pstpm2
summary(pstpm2(Surv(rectime,censrec==1)~hormon,data=brcancer2,weights=w))
summary(pstpm2(Surv(rectime,censrec==1)~hormon,data=brcancer2,weights=w,robust=TRUE))
##
## equal weights
brcancer2 <- transform(brcancer,w=4)
## unweighted
summary(fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer2,
      smooth.formula=~nsx(log(rectime),df=3,stata=TRUE)))
## weighted estimates
## stpm2
summary(stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer2,weights=w,
      smooth.formula=~nsx(log(rectime),df=3,stata=TRUE)))
summary(stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer2,weights=w,robust=TRUE,
      smooth.formula=~nsx(log(rectime),df=3,stata=TRUE)))
## pstpm2
summary(pstpm2(Surv(rectime,censrec==1)~hormon,data=brcancer2))
summary(pstpm2(Surv(rectime,censrec==1)~hormon,data=brcancer2,weights=w))
summary(pstpm2(Surv(rectime,censrec==1)~hormon,data=brcancer2,weights=w,robust=TRUE))

refresh
require(rstpm2)
brcancer2 <- transform(brcancer,w=ifelse(hormon==0,10,1))
##debug(rstpm2:::meat.stpm2)
summary(fit <- pstpm2(Surv(rectime,censrec==1)~hormon,data=brcancer2,weights=w,robust=TRUE,
      smooth.formula=~nsx(log(rectime),df=3,stata=TRUE)))
summary(fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer2,
      logH.formula=~nsx(log(rectime),df=3,stata=TRUE)))


## code for the SAS PROC ICPHREG examples
read.textConnection <- function(text, ...) {
    conn <-  textConnection(text)
    on.exit(close(conn))
    read.table(conn, ...)
}
hiv <- read.textConnection("0 16 0 0 0 1
15 26 0 0 0 1
12 26 0 0 0 1
17 26 0 0 0 1
13 26 0 0 0 1
0 24 0 0 1 0
6 26 0 1 1 0
0 15 0 1 1 0
14 26 0 1 1 0
12 26 0 1 1 0
13 26 0 1 0 1
12 26 0 1 1 0
12 26 0 1 1 0
0 18 0 1 0 1
0 14 0 1 0 1
0 17 0 1 1 0
0 15 0 1 1 0
3 26 1 0 0 1
4 26 1 0 0 1
1 11 1 0 0 1
13 19 1 0 0 1
0 6 1 0 0 1
0 11 1 1 0 0
6 26 1 1 0 0
0 6 1 1 0 0
2 12 1 1 0 0
1 17 1 1 1 0
0 14 1 1 0 0
0 25 1 1 0 1
2 11 1 1 0 0
0 14 1 1 0 0")
names(hiv) <- c("Left","Right","Stage","Dose","CdLow","CdHigh")
##hiv <- transform(hiv, Left=pmax(1e-5,Left))
hiv <- transform(hiv,Event = ifelse(Left==0,2,
                             ifelse(Right>=26,0,
                                    3)))
hiv2 <- transform(hiv,
                 Left = ifelse(Event==2,Right,
                        ifelse(Event==0,Left,
                               Left)),
                 Right = ifelse(Event==2,NA,
                         ifelse(Event==0,NA,
                                Right)))
library(rstpm2)
summary(stpm2(Surv(Left,Right,Event,type="interval")~Stage, data=hiv2, df=2))
summary(stpm2(Surv(Left,Right,Event,type="interval")~Stage, data=hiv2, df=1))
survreg(Surv(Left, Right, Event, type = "interval")~Stage, data=hiv2,dist="exponential")
library(rms)
psm(Surv(Left, Right, Event, type = "interval")~Stage, data=hiv2,dist="exponential")

##
library(rstpm2)
my.brcancer = brcancer
my.brcancer$left = my.brcancer$rectime
my.brcancer$right = ifelse(my.brcancer$censrec==1, my.brcancer$rectime, Inf)
test.stpm2.C = rstpm2::stpm2(Surv(left, right, censrec, type = "interval")~hormon,
                             data=my.brcancer,df=3)
## additive model
summary(fit <- stpm2(Surv(startTime,rectime,censrec==1)~hormon,data=brcancer2,
                     logH.formula=~nsx(rectime,df=3),
                     tvc.formula=~hormon:nsx(rectime,df=3,stata=TRUE)))

require(foreign)
require(rstpm2)
stmixed <- read.dta("http://fmwww.bc.edu/repec/bocode/s/stmixed_example2.dta")
system.time(r <- stpm2(Surv(stime,event)~treat,data=stmixed,cluster=stmixed$trial))
system.time(r <- stpm2(Surv(stime,event)~treat,data=stmixed,cluster=stmixed$trial,RandDist="LogN",
                       nodes=20))
summary(r)

require(mexhaz)
system.time(mix <-
                mexhaz(formula=Surv(stime,event)~treat, data=stmixed, base="exp.bs",degree=3,
                       random="trial", verbose=0))

## Frailty model
require(rstpm2)
require(frailtypack)
data(dataAdditive)
##debug(pstpm2)
system.time(mod2n <- pstpm2(Surv(t1,t2,event)~var1,
                           data=dataAdditive,
                           RandDist="LogN",
                           ##optimiser="NelderMead",
                           smooth.formula=~s(log(t2)),
                           sp.init=0.07723242,
                           adaptive=TRUE,
                           cluster=dataAdditive$group, nodes=10, trace=0))
summary(mod2n)

localargs <- mod2n@args
localargs$init <- mod2n@args$init*1.1
localargs$adaptive=TRUE
localargs$return_type <- "gradient"
.Call("model_output", localargs, package="rstpm2")
fdgrad <- function(obj,eps=1e-6) {
    args <- obj@args
    args$init <- args$init*1.1
    sapply(1:length(args$init), function(i) {
        args$return_type <- "objective"
        args$init[i] <- args$init[i]+eps
        f1 <- .Call("model_output", args, package="rstpm2")
        args$init[i] <- args$init[i]-2*eps
        f2 <- .Call("model_output", args, package="rstpm2")
        (f1-f2)/2/eps
    })
}
fdgrad(mod2n)
## OK for adaptive=FALSE

localargs <- mod2n@args
localargs$return_type <- "variances"
.Call("model_output", localargs, package="rstpm2")
localargs$return_type <- "modes"
.Call("model_output", localargs, package="rstpm2")

system.time(mod2nb <- stpm2(Surv(t1,t2,event)~var1,
                           data=dataAdditive,
                           RandDist="LogN",
                           logH.formula=~ns(log(t2),df=7),
                           cluster=dataAdditive$group, nodes=20))

system.time(mod2g <- pstpm2(Surv(t1,t2,event)~var1,
                           data=dataAdditive,
                           RandDist="Gamma",
                           smooth.formula=~s(log(t2)),
                           cluster=dataAdditive$group))


mod1 <- frailtyPenal(Surv(t1,t2,event)~cluster(group)+var1,data=dataAdditive,
                     n.knots=8,kappa=0.1,cross.validation=TRUE)
mod1n <- frailtyPenal(Surv(t1,t2,event)~cluster(group)+var1,data=dataAdditive,
                     n.knots=8,kappa=0.1,cross.validation=TRUE, RandDist="LogN")

system.time(mod2 <- stpm2(Surv(t1,t2,event)~var1, # Gamma
                          data=dataAdditive,
                          logH.formula=~ns(t2,df=7),
                          cluster=dataAdditive$group))

system.time(coxph1 <- coxph(Surv(t1,t2,event)~var1+frailty(group,distribution="gaussian"),
                          data=dataAdditive))
summary(coxph1)

system.time(mod2n <- stpm2(Surv(t1,t2,event)~var1,
                           data=dataAdditive,
                           RandDist="LogN",
                           optimiser="NelderMead",
                           logH.formula=~ns(log(t2),df=7),
                           cluster=dataAdditive$group, nodes=20))
system.time(mod2nb <- stpm2(Surv(t1,t2,event)~var1,
                           data=dataAdditive,
                           RandDist="LogN",
                           logH.formula=~ns(log(t2),df=7),
                           cluster=dataAdditive$group, nodes=20))

system.time(mod3 <- pstpm2(Surv(t1,t2,event)~var1,
                           data=dataAdditive,
                           RandDist="LogN",
                           criterion="BIC",
                           smooth.formula=~s(log(t2)),
                           cluster=dataAdditive$group, nodes=20))

system.time(mod3 <- coxph(Surv(t1,t2,event)~frailty(group,distribution="gamma")+var1,data=dataAdditive))
summary(mod2)
coef2 <- coef(summary(mod2))
theta <- exp(coef2[nrow(coef2),1])
se.logtheta <- coef2[nrow(coef2),2]
se.theta <- theta*se.logtheta
test.statistic <- 1/se.logtheta
pchisq(test.statistic,df=1,lower.tail=FALSE)/2

library(rstpm2)
library(ICE)
data(ICHemophiliac)
ICHemophiliac2 <- transform(as.data.frame(ICHemophiliac),event=3)
## fit1 <- stpm2(Surv(left,right,event,type="interval")~1,data=ICHemophiliac2,df=3)
fit1 <- pstpm2(Surv(left,right,event,type="interval")~1,data=ICHemophiliac2,
               smooth.formula=~s(left,k=7))
estimate <- ickde(ICHemophiliac, m=200, h=0.9)
plot(estimate, type="l", ylim=c(0,0.20))
tt <- seq(0,20,length=301)[-1]
lines(tt,predict(fit1,newdata=data.frame(left=tt),type="density"),col="blue")

## reg1 <- survreg(Surv(left,right,event,type="interval")~1,data=ICHemophiliac2)
## weibullShape <- 1/reg1$scale
## ## weibullScale <- exp(predict(reg1,type="lp"))
## weibullScale <- predict(reg1);
## tt <- seq(0,20,length=301)
## estimate <- ickde(ICHemophiliac, m=200, h=0.9)
## plot(estimate, type="l", ylim=c(0,0.15))
## lines(tt,dweibull(tt,weibullShape,weibullScale),lty=2)


library(rstpm2)
library(survival)
data(veteran)
## Re-define variables
veteran <- dplyr::mutate(veteran,
                         squamous = ifelse(celltype=="squamous",1,0),
                         smallcell = ifelse(celltype=="smallcell",1,0),
                         adeno = ifelse(celltype=="adeno",1,0),
                         large = ifelse(celltype=="large",1,0),
                         prior.ty = ifelse(prior==0,0,1),
                         trt = ifelse(trt==2,1,0),
                         high = ifelse(karno > 50,1,0))
lung<-subset(veteran, prior==0) ## patients with no prior therapy
## Why no optimal smoothing parameters?? divergence with version 1.3.3
pfit <-pstpm2(Surv(time,status==1) ~ adeno + smallcell + squamous,
                          smooth.formula = ~s(log(time)) + s(karno), data=lung, link.type="PO", trace = 1)

 

## two-dimensional smoothers
x1 <- x2 <- seq(0,1,length=11)
dat <- expand.grid(x1=x1,x2=x2)
dat$y <- rnorm(nrow(dat))
require(mgcv)
fit <- gam(y~s(x1,x2),data=dat)
fit$smooth


system.time(print(fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,type="probit")))
system.time(print(fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,type="probit",use.rcpp=FALSE)))

system.time(print(fit2 <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer, type="PO")))
system.time(print(fit2 <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer, type="PO",use.rcpp=FALSE)))

system.time(print(stpm2Gen(Surv(rectime,censrec==1)~hormon,data=brcancer)))
system.time(print(stpm2Gen(Surv(rectime,censrec==1)~hormon,data=brcancer, use.rcpp=FALSE)))
head(predict(fit,se.fit=TRUE))
head(predict(fit,type="haz",se.fit=TRUE))

brcancer <- brcancer[rep(1:nrow(brcancer),each=500),]
system.time(print(fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer))) # faster than Stata!
system.time(print(pfit <- pstpm2(Surv(rectime,censrec==1)~hormon,data=brcancer)))
plot(pfit,newdata=data.frame(hormon=0))

refresh
require(rstpm2)
data(brcancer)
system.time(print(fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,
                               tvc=list(hormon=3))))
system.time(print(pfit <- pstpm2(Surv(rectime,censrec==1)~1,data=brcancer,sp.init=c(0.0001,0.0001),
                                 tvc.formula=~s(log(rectime),by=hormon))))
print(pstpm2(Surv(rectime,censrec==1)~1,data=brcancer,init=coef(pfit)*100,
                                 tvc.formula=~s(log(rectime),by=hormon)))

summary(pfit)
plot(pfit,newdata=data.frame(hormon=0))
plot(pfit,newdata=data.frame(hormon=1),add=TRUE)
plot(pfit,newdata=data.frame(hormon=0),type="haz")
plot(pfit,newdata=data.frame(hormon=1),type="haz",add=TRUE)

pfit <- pstpm2(Surv(rectime/365,censrec==1)~1,data=brcancer) # OK
plot(pfit,newdata=data.frame(hormon=0))
system.time(print(pfit <- pstpm2(Surv(rectime/365,censrec==1)~1,data=brcancer,
               tvc.formula=~s(log(rectime/365),by=hormon))))
plot(pfit,newdata=data.frame(hormon=0)) # OK


times <- seq(500,2000,by=500)
meansurv1 <- t(sapply(times,function(time) predict(pfit,transform(brcancer,rectime=time,hormon=1),type="meansurv",se.fit=TRUE)))
meansurv0 <- t(sapply(times,function(time) predict(pfit,transform(brcancer,rectime=time,hormon=0),type="meansurv",se.fit=TRUE)))
matplot(times,meansurv1,type="l",lty=c(1,2,2),col=1)
matlines(times,meansurv0,type="l",lty=c(1,2,2),col=2)

meansurvdiff1 <- t(sapply(times,function(time) predict(pfit,transform(brcancer,rectime=time,hormon=0),type="meansurvdiff",var="hormon",se.fit=TRUE)))
matplot(times,meansurvdiff1,type="l",lty=c(1,2,2),col=1)



system.time(print(fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,control=list(parscale=100,reltol=1e-10),use.rcpp=FALSE)))

summary(fit)
system.time(print(fit2 <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,control=list(parscale=10000.0),reltol=1e-10,init=0.0001*coef(fit))))
summary(fit2)
plot(fit2,newdata=data.frame(hormon=1))

brcancerN <- brcancer[rep(1:nrow(brcancer),each=100),]
system.time(fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancerN,use.rcpp=FALSE,
                         control=list(parscale=0.1,reltol=1e-10)))
summary(fit)
system.time(fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancerN,use.rcpp=TRUE))
summary(fit)

###### penalised likelihood
## environment(pstpm2) <- environment(rstpm2::pstpm2)
## require(rstpm2)
try(detach("package:rstpm2",unload=TRUE))
## source("/home/MEB/marcle/src/R/rstpm2/R/pm2-3.R")

refresh
require(rstpm2)
data(brcancer)
brcancer$recyear <- brcancer$rectime/365
system.time(fit0 <- stpm2(Surv(recyear,censrec==1)~hormon,data=brcancer,df=5))
system.time(pfit0 <- pstpm2(Surv(recyear,censrec==1)~hormon,data=brcancer,sp.init=1))
system.time(pfit0.1 <- pstpm2(Surv(recyear,censrec==1)~hormon,data=brcancer,
                smooth.formula=~s(log(recyear),k=15),sp.init=10,alpha=2))
system.time(pfit1.1 <- pstpm2(Surv(recyear,censrec==1)~hormon,data=brcancer,
                smooth.formula=~s(log(recyear)),sp.init=10,criterion="BIC"))
system.time(pfit2 <- pstpm2(Surv(recyear,censrec==1)~hormon,data=brcancer,
                smooth.formula=~s(recyear),sp.init=10))


plot(pfit0,newdata=data.frame(hormon=1),line.col="red",type="hazard")
plot(pfit0.1,newdata=data.frame(hormon=1),line.col="blue",add=TRUE,type="hazard")
plot(pfit1.1,newdata=data.frame(hormon=1),line.col="orange",add=TRUE,type="hazard")
plot(fit0,newdata=data.frame(hormon=1),line.col="green",type="hazard",add=TRUE)
plot(pfit2,newdata=data.frame(hormon=1),line.col="black",type="hazard",add=TRUE)

plot(pfit0,newdata=data.frame(hormon=1),line.col="red")
plot(pfit0.1,newdata=data.frame(hormon=1),line.col="blue",add=TRUE)
plot(pfit1.1,newdata=data.frame(hormon=1),line.col="pink",add=TRUE)
plot(fit0,newdata=data.frame(hormon=1),line.col="green",add=TRUE)
plot(pfit2,newdata=data.frame(hormon=1),line.col="black",add=TRUE)



system.time(pfit0.check <- pstpm2(Surv(recyear,censrec==1)~hormon,data=brcancer, sp=pfit0@sp, use.rcpp=FALSE))
system.time(pfit0.check2 <- pstpm2(Surv(recyear,censrec==1)~hormon,data=brcancer, sp=pfit0@sp))
summary(pfit0)
summary(pfit0.check)
summary(pfit0.check2)

system.time(pfit0 <- pstpm2(Surv(recyear,censrec==1)~hormon,data=brcancer,
                logH.formula=~s(log(recyear),k=30),sp.init=1))
system.time(pfit0.1 <- pstpm2(Surv(recyear,censrec==1)~hormon,data=brcancer,
                logH.formula=~s(log(recyear),k=30),sp.init=1,criterion="BIC"))
plot(pfit0,newdata=data.frame(hormon=1),line.col="red",type="hazard")
plot(pfit0.1,newdata=data.frame(hormon=1),line.col="blue",add=TRUE,type="hazard")

system.time(pfit1 <- pstpm2(Surv(recyear,censrec==1)~hormon,data=brcancer,
                logH.formula=~s(log(recyear),k=30),sp=10,pen="h",
                            smoother.parameters=list("log(recyear)"=list(var="recyear",
                                                         inverse=exp,
                                                         transform=log))))
plot(pfit1,newdata=data.frame(hormon=1),line.col="green",add=TRUE,type="hazard")


system.time(pfit1 <- pstpm2(Surv(recyear,censrec==1)~hormon,data=brcancer,
                logH.formula=~s(log(recyear),k=30),sp=10,criterion="GCV"))
system.time(pfit2 <- pstpm2(Surv(recyear,censrec==1)~hormon,data=brcancer,
                logH.formula=~s(log(recyear),k=20),sp=10,criterion="BIC"))
plot(pfit1,newdata=data.frame(hormon=1),type="hazard",ylim=c(0,0.25))
plot(pfit2,newdata=data.frame(hormon=1),add=TRUE,line.col="blue",type="hazard")

system.time(pfit1 <- pstpm2(Surv(recyear,censrec==1)~hormon,data=brcancer,
                logH.formula=~s(log(recyear),k=30),sp=1,criterion="GCV"))

system.time(pfit2 <- pstpm2(Surv(recyear,censrec==1)~hormon,data=brcancer,
                logH.formula=~s(recyear,k=20),sp=1,use.rcpp=F,penalty="h"))
system.time(pfit2 <- pstpm2(Surv(recyear,censrec==1)~hormon,data=brcancer,
                logH.formula=~s(recyear,k=20),sp=1,penalty="h"))
plot(pfit2,newdata=data.frame(hormon=1),type="hazard")

system.time(pfit2 <- pstpm2(Surv(recyear,censrec==1)~hormon,data=brcancer,
                logH.formula=~s(recyear,k=30),sp=1,use.rcp=FALSE))

plot(pfit1,newdata=data.frame(hormon=1),type="hazard")
plot(pfit2,newdata=data.frame(hormon=1),type="hazard")

system.time(pfit2.0 <- pstpm2(Surv(recyear,censrec==1)~hormon,data=brcancer,
                logH.formula=~s(recyear,k=30),sp=0.055,penalty="h",cr="GCV"))
system.time(pfit2 <- pstpm2(Surv(recyear,censrec==1)~hormon,data=brcancer,
                logH.formula=~s(recyear,k=30),sp=0.055,use.rcp=FALSE,penalty="h"))
rstpm2:::gcv(pfit2)
plot(pfit2,newdata=data.frame(hormon=1),line.col="red",add=TRUE,type="hazard")
plot(pfit2.0,newdata=data.frame(hormon=1),line.col="green",add=TRUE,type="hazard")

require(frailtypack)
fpack1 <- frailtyPenal(Surv(recyear,censrec==1)~hormon, data=brcancer, cross.validation=TRUE, n.knots=10, kappa1=0.1)
plot(fpack1)

system.time(pfit1 <- pstpm2(Surv(recyear,censrec==1)~hormon+x3,data=brcancer,
                logH.formula=~s(log(recyear),k=20)+s(x3),sp=c(0.1,0.1)))
system.time(pfit1 <- pstpm2(Surv(recyear,censrec==1)~hormon+x3,data=brcancer,
                logH.formula=~s(log(recyear),k=20)+s(x3)))
system.time(pfit1 <- pstpm2(Surv(recyear,censrec==1)~hormon,data=brcancer,
                logH.formula=~s(log(recyear),k=20)+s(x3)))
plot(pfit1,newdata=data.frame(hormon=1,x3=20))
plot(pfit1,newdata=data.frame(hormon=0,x3=20),type="hazard")
plot(pfit1,newdata=data.frame(hormon=1,x3=20),type="hazard",add=TRUE,line.col="blue",lty=1)

summary(pfit1)
brcancerN <- brcancer[rep(1:nrow(brcancer),each=100),]
system.time(pfit1 <- pstpm2(Surv(recyear,censrec==1)~hormon,data=brcancerN,
                logH.formula=~s(log(recyear),k=20)))
plot(pfit1,newdata=data.frame(hormon=1))
pfit1@gam$sp
par(mfrow=c(2,2))
plot(pfit1,newdata=data.frame(hormon=1))
summary(pfit1@gam)$edf
rstpm2:::gcv(pfit1)
rstpm2:::gcvc(pfit1,nn)
sps <- as.list(10^(seq(-4,2,by=0.5)))
system.time(pfit2 <- pstpm2(Surv(recyear,censrec==1)~hormon,data=brcancer,
                logH.formula=~s(log(recyear),k=20), sp=sps))
gcvs <- lapply(pfit2,rstpm2:::gcv)
plot(sps,unlist(gcvs),type="l",log="x")
plot(sapply(gcvs,attr,"negll"),sapply(gcvs,attr,"trace"),type="l",asp=1)
plot(sapply(gcvs,attr,"trace"),sapply(gcvs,attr,"negll"),type="l",asp=1)
plot(sps,sapply(pfit2,rstpm2:::aicc,nn=nn),type="l",log="x")
plot(sps,sapply(pfit2,rstpm2:::bic,nn=nn),type="l",log="x")
##gcvc
brcancer$recyear <- brcancer$rectime/365
sps <- 10^(seq(-4,2,by=0.5))
gcvcs <- sapply(sps, function(sp) {
gcvc(pstpm2(Surv(recyear,censrec==1)~hormon,data=brcancer,
                      logH.formula=~s(recyear,k=30), sp=sp),length(brcancer$recyear))
})
plot(sps,gcvcs,type="l",log="x")
###bic
brcancer$recyear <- brcancer$rectime/365
sps <- 10^(seq(-4,2,by=0.5))
gcvcs <- sapply(sps, function(sp) {
  bic(pstpm2(Surv(recyear,censrec==1)~hormon,data=brcancer,
              logH.formula=~s(recyear,k=30), sp=sp),length(brcancer$recyear))
})
plot(sps,gcvcs,type="l",log="x")
###aicc
brcancer$recyear <- brcancer$rectime/365
sps <- 10^(seq(-4,2,by=0.5))
gcvcs <- sapply(sps, function(sp) {
  aicc(pstpm2(Surv(recyear,censrec==1)~hormon,data=brcancer,
             logH.formula=~s(recyear,k=30), sp=sp),length(brcancer$recyear))
})
plot(sps,gcvcs,type="l",log="x")
#########################


### penalty functions
require(mgcv)
require(gaussquad)

## Outline:
## get w, lambda, X0, X1, X2, X3
## calculate s0, s1, s2, s3
## calculate h2 and pfun=integrate(h2^2,t)
## calculate dh2sq.dbeta and dpfun=integrate(dh2sq.dbeta,t)
##
## calculate w, lambda, X0, X1, X2, X3
derivativeDesign <- 
function (functn, lower = -1, upper = 1, rule = NULL,
    ...) 
{
    pred <- if (length(list(...)) && length(formals(functn)) > 
              1) 
        function(x) functn(x, ...)
    else functn
    if (is.null(rule))
        rule <-    ## gaussquad::legendre.quadrature.rules(20)[[20]]
        data.frame(x = c(0.993128599185095, 0.963971927277914, 0.912234428251326, 
                       0.839116971822219, 0.746331906460151, 0.636053680726515, 0.510867001950827, 
                       0.37370608871542, 0.227785851141646, 0.0765265211334977, -0.0765265211334974, 
                       -0.227785851141645, -0.373706088715418, -0.510867001950827, -0.636053680726516, 
                       -0.746331906460151, -0.839116971822219, -0.912234428251326, -0.963971927277913, 
                       -0.993128599185094),
                   w = c(0.0176140071391522, 0.040601429800387, 
                       0.0626720483341092, 0.0832767415767053, 0.101930119817241, 0.11819453196152, 
                       0.131688638449176, 0.14209610931838, 0.149172986472603, 0.152753387130726, 
                       0.152753387130726, 0.149172986472603, 0.142096109318381, 0.131688638449175, 
                       0.11819453196152, 0.10193011981724, 0.0832767415767068, 0.0626720483341075, 
                       0.0406014298003876, 0.0176140071391522))
    lambda <- (upper - lower)/(2)
    mu <- (lower + upper)/(2)
    x <- lambda * rule$x + mu
    w <- rule$w
    eps <- .Machine$double.eps^(1/8)
    X0 <- pred(x)
    X1 <- (-pred(x+2*eps)+8*pred(x+eps)-8*pred(x-eps)+pred(x-2*eps))/12/eps
    X2 <- (-pred(x+2*eps)/12+4/3*pred(x+eps)-5/2*pred(x)+4/3*pred(x-eps)-pred(x-2*eps)/12)/eps/eps
    X3 <- (-pred(x+3*eps)/8+pred(x+2*eps)-13/8*pred(x+eps)+
           13/8*pred(x-eps)-pred(x-2*eps)+pred(x-3*eps)/8)/eps/eps/eps
    return(list(x=x,w=w,lambda=lambda,X0=X0,X1=X1,X2=X2,X3=X3))
}
hpfun <- function(beta,design) {
    lapply(design,function(obj) {
        s0 <- as.vector(obj$X0 %*% beta)
        s1 <- as.vector(obj$X1 %*% beta)
        s2 <- as.vector(obj$X2 %*% beta)
        s3 <- as.vector(obj$X3 %*% beta)
        h2 <- (s3+3*s1*s2+s1^3)*exp(s0)
        obj$lambda*sum(obj$w*h2^2)
    })
}
hdpfun <- function(beta,design) {
    lapply(design, function(obj) {
        s0 <- as.vector(obj$X0 %*% beta)
        s1 <- as.vector(obj$X1 %*% beta)
        s2 <- as.vector(obj$X2 %*% beta)
        s3 <- as.vector(obj$X3 %*% beta)
        h2 <- (s3+3*s1*s2+s1^3)*exp(s0)
        dh2sq.dbeta <- 2*h2*(exp(s0)*(obj$X3+3*(obj$X1*s2+obj$X2*s1)+3*s1^2*obj$X1)+h2*obj$X0)
        obj$lambda*colSums(obj$w*dh2sq.dbeta)
    })
}
smootherDesign <- function(gamobj,data) {
    d <- data[1,,drop=FALSE] ## how to get mean prediction values, particularly for factors?
    makepred <- function(var) {
        function(value) {
            d <- d[rep(1,length(value)),]
            d[[var]] <- value
            predict(gamobj,newdata=d,type="lpmatrix")
        }
    }
    lapply(gamobj$smooth, function(smoother) {
        var <- smoother$term
        pred <- makepred(var)
        derivativeDesign(pred,lower=min(data[[var]]),upper=max(data[[var]]))
    })
}
## example data
d <- within(data.frame(x=seq(0,1,length=301)), {
    mu <- exp(x)
    y <- rnorm(301,mu,0.01)
})
fit <- gam(y~s(x),data=d,family=gaussian(link="log"))
beta <- coef(fit)
design <- smootherDesign(fit,d)
hpfun(beta,design)
hdpfun(beta,design)

    

## Testing...
require(mgcv)
d <- within(data.frame(x=seq(0,2*pi,length=301)), {
    mu <- sin(x)
    dmu <- cos(x)
    y <- rnorm(301,mu,0.001)
})
fit <- gam(y~s(x),data=d)
mat <- predict(fit,type="lpmatrix")
with(d,plot(x,y))
with(d,lines(x,mu,lwd=2))
with(d,lines(x,predict(fit),col="blue",lwd=2))
par(mfrow=c(3,2))
pred <- function(eps,obj=fit,data=d,var="x") {
    nd <- d
    nd[[var]] <- nd[[var]]+eps
    predict(obj,newdata=nd,type="lpmatrix")
}
## First derivative
eps <- .Machine$double.eps^(1/8)
matD <- (pred(eps) - pred(-eps)) / 2 / eps
with(d,plot(x,dmu,lwd=2,type="l"))
with(d,lines(x,matD %*% coef(fit),col="blue",lwd=2))
##
## 1/12 	−2/3 	0 	2/3 	−1/12
eps <- .Machine$double.eps^(1/8)
matD <- (-pred(2*eps)+8*pred(eps)-8*pred(-eps)+pred(-2*eps))/12/eps
with(d,plot(x,dmu,lwd=2,type="l"))
with(d,lines(x,matD %*% coef(fit),col="blue",lwd=2))
##
## Second derivative
eps <- .Machine$double.eps^(1/8)
matD2 <- (pred(eps)-2*pred(0)+pred(-eps))/eps/eps
with(d,plot(x,-mu,lwd=2,type="l"))
with(d,lines(x,matD2 %*% coef(fit),col="blue",lwd=2))
##
## −1/12 	4/3 	−5/2 	4/3 	−1/12 	
eps <- .Machine$double.eps^(1/8)
matD2 <- (-pred(2*eps)/12+4/3*pred(eps)-5/2*pred(0)+4/3*pred(-eps)-pred(-2*eps)/12)/eps/eps
with(d,plot(x,-mu,lwd=2,type="l"))
with(d,lines(x,matD2 %*% coef(fit),col="blue",lwd=2))
##
## Third derivatives
eps <- .Machine$double.eps^(1/8)
matD3 <- (pred(2*eps)-
          2*pred(eps)+
          2*pred(-eps)-
          pred(-2*eps))/2/eps/eps/eps
with(d,plot(x,-dmu,lwd=2,type="l"))
with(d,lines(x,matD3 %*% coef(fit),col="blue",lwd=2))
##
## 1/8 	−1 	13/8 	0 	−13/8 	1 	−1/8
eps <- .Machine$double.eps^(1/8)
matD3 <- (-pred(3*eps)/8+pred(2*eps)-13/8*pred(eps)+
          13/8*pred(-eps)-pred(-2*eps)+pred(-3*eps)/8)/eps/eps/eps
with(d,plot(x,-dmu,lwd=2,type="l"))
with(d,lines(x,matD3 %*% coef(fit),col="blue",lwd=2))


## (browse-url "http://en.wikipedia.org/wiki/Finite_difference_coefficients")



require(mgcv)
data <- data.frame(x=1:10,y=1:10)
fit <- gam(y~s(x,k=5,bs="ps"),data=data)

round(cbind(1,(spline.des(knots=fit$smooth[[1]]$knots,x=data$x)$design %*% 
  qr.Q(attr(fit$smooth[[1]],"qrc"),complete=TRUE))[,-1]) -
           predict(fit,type="lpmatrix"),1e-10)

round(cbind(1,(spline.des(knots=fit$smooth[[1]]$knots,x=5:6)$design %*% 
                 qr.Q(attr(fit$smooth[[1]],"qrc"),complete=TRUE))[,-1]) -
        predict(fit,newdata=data.frame(x=5:6),type="lpmatrix"),1e-10)

cbind(0,(spline.des(knots=fit$smooth[[1]]$knots,x=data$x,derivs=rep(1,nrow(data)))$design %*% 
        qr.Q(attr(fit$smooth[[1]],"qrc"),complete=TRUE))[,-1]) -
  (predict(fit,newdata=transform(data,x=x+1e-5),type="lpmatrix")-
     (predict(fit,newdata=transform(data,x=x-1e-5),type="lpmatrix")))/2e-5



#######Optimal fitting#######
###GCV,AICC,BIC or GCVC to choose smoothing parameters###
opt.val<-function(pstpm2.fit,nn){
  like<-pstpm2.fit@like
  Hl<-numDeriv::hessian(like,coef(pstpm2.fit))
  Hinv<-vcov(pstpm2.fit)
  trace<-sum(diag(Hinv%*%Hl))
  loglike<-(like(coef(pstpm2.fit)))/nn
  gcv<-(trace-loglike)/nn
  aicc<-(-2*loglike+2*trace*nn/(nn-trace-1))/nn
  bic<-(-2*loglike+trace*log(nn))/nn
  gcvc<-(-2*loglike-2*nn*log(1-trace/nn))/nn
  out<-c(loglike,gcv,aicc,bic,gcvc)
  return(out)
}
###############################
###############################
# setClass("opt.fit", representation(
#                                    num.ind = "numeric",
#                                    cr = "numeric",
#                                    tops = "data.frame",
#                                    sp.opt = "numeric",
#                                    fun.min = "numeric"
# ),
#          contains="pstpm2")
# #########################
# opt.fit<-function(formula,data,logH.formula,sp.low,sp.upp,num.sp,timeVar = NULL){
#   ###number of individual
#   num.ind <- nrow(data)
#   #####Censoring rate####
#   ## set up the data
#   ## ensure that data is a data frame
#   data <- get_all_vars(formula, data)
#   #   ## parse the function call
#   #   Call <- match.call()
#   #   mf <- match.call(expand.dots = FALSE)
#   #   m <- match(c("formula", "data", "subset", "contrasts", "weights"),
#   #              names(mf), 0L)
#   #   mf <- mf[c(1L, m)]
#   stopifnot(length(lhs(formula))>=2)
#   eventExpr <- lhs(formula)[[length(lhs(formula))]]
#   delayed <- length(lhs(formula))==4
#   timeExpr <- lhs(formula)[[if (delayed) 3 else 2]] # expression
#   if (is.null(timeVar))
#     timeVar <- all.vars(timeExpr)
#   time <- eval(timeExpr, data)
#   if (delayed) {
#     time0Expr <- lhs(formula)[[2]]
#     time0 <- eval(time0Expr, data)
#   }
#   event <- eval(eventExpr,data)
#   cr <- sum(event > min(event))/num.ind
#   #   
#   #   cr=table(lhs(formula)[[if (delayed) 4 else 3]][2])/nn
#   ##nn<-length(brcancer$recyear)
#   # system.time(pfit1 <- pstpm2(Surv(recyear,censrec==1)~hormon,data=brcancer,
#   #                             logH.formula=~s(recyear,k=30), sp=1e-1))
#   # plot(pfit1,newdata=data.frame(hormon=1))
#   
#   #sps <- 10^(seq(-4,4,by=0.5))
#   #   sp.low=10^-4
#   #   sp.upp=4000
#   #   num.sp=30
#   sps <- 10^(seq(log10(sp.low),log10(sp.upp),length=num.sp))
#   optvals <- sapply(sps, function(sp) {
#     opt.val(pstpm2(formula,data,logH.formula=NULL, sp=sp),num.ind)
#   })
#   tops<-t(optvals)
#   colnames(tops) <- c("loglike","gcv","aicc","bic","gcvc")
#   rownames(tops) <- rownames(tops, do.NULL = FALSE, prefix = "Obs.")
#   # tops<-as.data.frame(tops)
#   tops<-as.data.frame(tops)
#   ####Plot#########
#   #par(mfrow=c(1,2))
#   ###to choose optimal smoothing parameter ###
#   ind.min <- sapply(2:5,function(x) order(tops[,x])[1])
#   sp.opt <- sps[ind.min]
#   obj<-pstpm2(formula,data,logH.formula=NULL, sp=sp.opt[1])
#   fun.min <- sapply(2:5,function(x) min(tops[,x]))
#   # if(ind.min[1]==1)
#   # stop("Hit left boundary, make sp.low smaller.")
#   # if(ind.min[1]==num.sp)
#   # stop("Hit right boundary, make sp.upp bigger.")
#   #   with(tops,matplot(sps,tops[,-1],type="l",col=1:4,lty=1:4,xlab="x",ylab="y"))
#   #   points(sp.opt,fun.min,pch=4,lwd=2,cex=1.2)
#   #   lines(sp.opt,fun.min,err=-1,col=1:4,lty=1:4)
#   
#   ###Estimate final model with optimal value of sp###
#   
#   #   
#   #   summary(pfit.obj)
#   #########################################
#   out <- as(obj,"opt.fit")
#   out <- new("opt.fit",
#              coef = pstpm2@coef,
#              fullcoef = pstpm2@fullcoef,
#              vcov = pstpm2@vcov,
#              min = pstpm2@min,
#              details = pstpm2@details,
#              minuslogl = pstpm2@minuslogl,
#              method = pstpm2@method,
#              data = data,
#              formula = pstpm2@formula,
#              optimizer = "optim",
#              xlevels = .getXlevels(pstpm2@terms,pstpm2@model.frame),
#              ##contrasts = attr(X, "contrasts"),
#              contrasts = NULL, # wrong!
#              logli = pstpm2@logli,
#              ##weights = weights,
#              Call = pstpm2@Call,
#              terms = pstpm2@terms,
#              model.frame = pstpm2@model.frame,
#              gam = pstpm2@gam,
#              timeVar = pstpm2@timeVar,
#              timeExpr = pstpm2@timeExpr,
#              like = pstpm2@like,
#              negll<-pstpm2@negll,
#              call.formula = pstpm2@call.formula,
#              x = pstpm2@x,
#              xd = pstpm2@xd,
#              termsd = pstpm2@termsd, # wrong!
#              y = pstpm2@y,
#              num.ind = num.ind,
#              cr = cr,
#              tops = tops,
#              sp.opt = sp.opt,
#              fun.min = fun.min)
#   
#   return(out)
# }


#####load data####
load("brcancer.rda")
data(brcancer)
brcancer$recyear <- brcancer$rectime/365
####model fit###
opt.fit(Surv(recyear,censrec==1)~hormon,data=brcancer,
        logH.formula=~s(recyear), sp.low=10^-4,sp.upp=4000,
        num.sp=30,timeVar = NULL)
# ###methods for Plot ###
# setMethod(
#   f= "plot",
#   signature(x="opt.fit", y="missing"),
#   definition=function (x,y,...){
#     matplot(x@sps,x@tops[,-1],type="l",col=1:4,lty=1:4,xlab="",ylab="")
#     points(x@sp.opt,x@fun.min,pch=4,lwd=2,cex=1.2)
#     lines(x@sp.opt,x@fun.min,err=-1,col=1:4,lty=1:4)
#   }
# )
# ####methods for print####
# setMethod ("print",signature(x="opt.fit", y="missing"),
#            function(x,...){
#              cat("*** Class opt.fit, method Print *** \n")
#              cat("* Optimal SP ="); print (x@sp.opt)
#              cat("* GCV = \n"); print (x@fun.min[1])
#              cat("******* End Print (opt.fit) ******* \n")
#            }
# )
##########################

aplot <- 
function (x, y, ...) 
{
    .local <- function (x, y, newdata, type = "surv", xlab = NULL, 
        ylab = NULL, line.col = 1, ci.col = "grey", lty = par("lty"), 
        add = FALSE, ci = !add, rug = !add, var = NULL, ...) 
    {
        browser()
        y <- predict(x, newdata, type = type, var = var, grid = TRUE, 
            se.fit = TRUE)
        if (is.null(xlab)) 
            xlab <- deparse(x@timeExpr)
        if (is.null(ylab)) 
            ylab <- switch(type, hr = "Hazard ratio", hazard = "Hazard", 
                surv = "Survival", sdiff = "Survival difference", 
                hdiff = "Hazard difference", cumhaz = "Cumulative hazard")
        xx <- attr(y, "newdata")
        xx <- eval(x@timeExpr, xx)
        if (!add) 
            matplot(xx, y, type = "n", xlab = xlab, ylab = ylab, 
                ...)
        if (ci) 
            polygon(c(xx, rev(xx)), c(y[, 2], rev(y[, 3])), col = ci.col, 
                border = ci.col)
        lines(xx, y[, 1], col = line.col, lty = lty)
        if (rug) {
            Y <- x@y
            eventTimes <- Y[Y[, ncol(Y)] == 1, ncol(Y) - 1]
            rug(eventTimes, col = line.col)
        }
        return(invisible(y))
    }
    .local(x, y, ...)
}
aplot(fit,newdata=data.frame(hormon=1))

apredict <- function (object, ...) 
{
    .local <- function (object, newdata = NULL, type = c("surv", 
        "cumhaz", "hazard", "hr", "sdiff", "hdiff", "loghazard", 
        "link"), grid = FALSE, seqLength = 300, se.fit = FALSE, 
        link = NULL, exposed = incrVar(var), var, ...) 
    {
        local <- function(object, newdata = NULL, type = "surv", 
            exposed) {
            ## browser()
            tt <- object@terms
            if (is.null(newdata)) {
                X <- object@x
                XD <- object@xd
                y <- object@y
                time <- as.vector(y[, ncol(y) - 1])
            }
            else {
                lpfunc <- function(delta, fit, data, var) {
                  data[[var]] <- data[[var]] + delta
                  lpmatrix.lm(fit, data)
                }
                X <- lpmatrix.lm(object@lm, newdata)
                XD <- grad(lpfunc, 0, object@lm, newdata, object@timeVar)
                XD <- matrix(XD, nrow = nrow(X))
                if (type %in% c("hazard", "hr", "sdiff", "hdiff", 
                  "loghazard")) {
                  time <- eval(object@timeExpr, newdata)
                }
                if (object@delayed) {
                  newdata0 <- newdata
                  newdata0[[object@timeVar]] <- newdata[[object@time0Var]]
                  X0 <- lpmatrix.lm(object@lm, newdata0)
                }
                if (type %in% c("hr", "sdiff", "hdiff")) {
                  if (missing(exposed)) 
                    stop("exposed needs to be specified for type in ('hr','sdiff','hdiff')")
                  newdata2 <- exposed(newdata)
                  X2 <- lpmatrix.lm(object@lm, newdata2)
                  XD2 <- grad(lpfunc, 0, object@lm, newdata2, 
                    object@timeVar)
                  XD2 <- matrix(XD, nrow = nrow(X))
                }
            }
            beta <- coef(object)
            cumHaz = as.vector(exp(X %*% beta))
            Sigma = vcov(object)
            if (type == "link") {
                return(as.vector(X %*% beta))
            }
            if (type == "cumhaz") {
                if (object@delayed) 
                  return(cumHaz - as.vector(X0 %*% beta))
                else return(cumHaz)
            }
            if (type == "surv") {
                return(exp(-cumHaz))
            }
            if (type == "sdiff") 
                return(as.vector(exp(-exp(X2 %*% beta))) - exp(-cumHaz))
            if (type == "hazard") {
                return(as.vector(XD %*% beta) * cumHaz)
            }
            if (type == "loghazard") {
                return(as.vector(log(XD %*% beta)) + log(cumHaz))
            }
            if (type == "hdiff") {
                return(as.vector((XD2 %*% beta) * exp(X2 %*% beta) - (XD %*% 
                  beta)/time * cumHaz))
            }
            if (type == "hr") {
                cumHazRatio = exp((X2 - X) %*% beta)
                return(as.vector((XD2 %*% beta)/(XD %*% beta) * cumHazRatio))
            }
        }
        type <- match.arg(type)
        if (is.null(newdata) && type %in% c("hr", "sdiff", "hdiff")) 
            stop("Prediction using type in ('hr','sdiff','hdiff') requires newdata to be specified.")
        if (grid) {
            Y <- object@y
            event <- Y[, ncol(Y)] == 1
            time <- object@data[[object@timeVar]]
            eventTimes <- time[event]
            X <- seq(min(eventTimes), max(eventTimes), length = seqLength)[-1]
            data.x <- data.frame(X)
            names(data.x) <- object@timeVar
            newdata <- merge(newdata, data.x)
        }
        pred <- if (!se.fit) {
            local(object, newdata, type = type, exposed = exposed, 
                ...)
        }
        else {
            if (is.null(link)) 
                link <- switch(type, surv = "cloglog", cumhaz = "log", 
                  hazard = "log", hr = "log", sdiff = "I", hdiff = "I", 
                  loghazard = "I", link = "I")
            predictnl(object, local, link = link, newdata = newdata, 
                type = type, exposed = exposed, ...)
        }
        attr(pred, "newdata") <- newdata
        return(pred)
    }
    .local(object, ...)
}
environment(apredict) <- environment(stpm2)
dim(apredict(fit,newdata=data.frame(hormon=1),grid=T)) # n=300 or 299??
apredict(fit,newdata=data.frame(hormon=1),grid=T)
dim(apredict(fit,newdata=data.frame(hormon=1),grid=T,se.fit=T)) # n=300 or 299??
apredict(fit,newdata=data.frame(hormon=1),grid=T,se.fit=T)

debug(rstpm2:::numDeltaMethod)

try(suppressWarnings(detach("package:rstpm2",unload=TRUE)),silent=TRUE)
require(rstpm2)
data(brcancer)
system.time(fit2 <- stpm2(Surv(rectime/365,censrec==1)~hormon,data=brcancer,df=5))
system.time(fit3 <- pstpm2(Surv(rectime/365,censrec==1)~hormon,data=brcancer,use.gr=F))
plot(fit3,newdata=data.frame(hormon=0),type="hazard")
plot(fit2,newdata=data.frame(hormon=0),type="hazard",add=TRUE,ci=FALSE,rug=FALSE,
     line.col=2)

## penalised likelihood
brcancer$recyear <- brcancer$rectime/365
system.time(pfit1 <- pstpm2(Surv(recyear,censrec==1)~hormon,data=brcancer,
                            logH.formula=~s(log(recyear),k=30), sp=1e-1))
system.time(fit1 <- stpm2(Surv(recyear,censrec==1)~hormon,data=brcancer,
                            logH.formula=~ns(log(recyear),df=4)))
plot(pfit1,newdata=data.frame(hormon=1))
plot(fit1,newdata=data.frame(hormon=1),lty=2,add=TRUE,ci=F)
rstpm2:::gcv(pfit1)
sps <- 10^(seq(-4,2,by=0.5))
gcvs <- sapply(sps, function(sp) {
  rstpm2:::gcv(pstpm2(Surv(recyear,censrec==1)~hormon,data=brcancer,
        logH.formula=~s(recyear,k=30), sp=sp))
})
plot(sps,gcvs,type="l",log="x")
##
system.time(fit <- rstpm2:::stpm2Old(Surv(rectime,censrec==1)~hormon,df=5,data=brcancer))
system.time(fit2 <- stpm2(Surv(rectime,censrec==1)~hormon,df=5,data=brcancer))
system.time(fit3 <- pstpm2(Surv(rectime,censrec==1)~hormon,data=brcancer))
##
plot(fit3,newdata=data.frame(hormon=0),type="hazard")
plot(fit2,newdata=data.frame(hormon=0),type="hazard",add=TRUE,line.col=2,ci=FALSE)
##
system.time(fit <- stpm2(Surv(rectime/365,censrec==1)~hormon,df=5,data=brcancer))
system.time(fit2 <- rstpm2:::stpm2Old(Surv(rectime/365,censrec==1)~hormon,df=5,data=brcancer))
##
system.time(fit3 <- pstpm2(Surv(rectime/365,censrec==1)~hormon,data=brcancer))
plot(fit3,newdata=data.frame(hormon=0),type="hazard")
##
plot(fit2,newdata=data.frame(hormon=0),type="hazard",add=TRUE,line.col=2,ci=FALSE)
##
plot(fit.tvc,newdata=data.frame(hormon=1),type="hr",var="hormon")
##
summary(fit.tvc <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,df=3,
                         tvc=list(hormon=3)))
anova(fit,fit.tvc) # compare with and without tvc
summary(fit.tvc <- stpm2Old(Surv(rectime,censrec==1)~hormon,data=brcancer,df=3,
                                      tvc=list(hormon=3)))
anova(fit,fit.tvc) # compare with and without tvc
##
plot(fit.tvc,newdata=data.frame(hormon=0),type="hr",var="hormon")
                                        # no lines method: use add=TRUE
plot(fit.tvc,newdata=data.frame(hormon=1),type="hr",var="hormon",
     add=TRUE,ci=FALSE,line.col=2)
##
## plain: identical results (good)
stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer)
stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,
                         logH.formula=~ns(log(rectime),3))
rstpm2:::stpm2Old(Surv(rectime,censrec==1)~hormon,data=brcancer)
## cure: identical (requires bhazard to be sensible)
rate0 <- 10^(-5+brcancer$x1/100)
(fit1 <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,df=2,cure=T,bhazard=rate0))
(fit2 <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,
                         logH.formula=~nsx(log(rectime),df=2,cure=T,log=T),bhazard=rate0))
(fit3 <- rstpm2:::stpm2Old(Surv(rectime,censrec==1)~hormon,data=brcancer,cure=T,df=2,bhazard=rate0))
(fit4 <- rstpm2:::stpm2Old(Surv(rectime,censrec==1)~hormon,data=brcancer,bhazard=rate0,
                           logH.formula=~nsx(log(rectime),2,cure=T)))


##### examples #####
require(foreign)
if (FALSE) { # testing in open code
  install.packages("bbmle", repos="http://R-Forge.R-project.org")
  require(bbmle)
  brcancer=read.dta("brcancer.dta")
  brcancer=transform(brcancer,rate0=10^(-5+x1/100))
}
try(suppressWarnings(detach("package:bbmle",unload=TRUE)),silent=TRUE)

try(suppressWarnings(detach("package:rstpm2",unload=TRUE)),silent=TRUE)
## require(rstpm2)
data(brcancer)
summary(fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,
                     logH.formula=~nsx(log(rectime),df=3,stata=TRUE)))

brcancer <- transform(brcancer,w=10)
summary(fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,
                     weights=w,robust=TRUE,
                     logH.formula=~nsx(log(rectime),df=3,stata=TRUE)))


## sandwich variance estimator (from the sandwich package)

coeftest.stpm2 <- 
function (x, vcov. = NULL, df = NULL, ...) 
{
    est <- coef(x)
    if (is.null(vcov.)) 
        se <- vcov(x)
    else {
        if (is.function(vcov.)) 
            se <- vcov.(x)
        else se <- vcov.
    }
    se <- sqrt(diag(se))
    if (!is.null(names(est)) && !is.null(names(se))) {
        anames <- names(est)[names(est) %in% names(se)]
        est <- est[anames]
        se <- se[anames]
    }
    tval <- as.vector(est)/se
    pval <- 2 * pnorm(abs(tval), lower.tail = FALSE)
    cnames <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    mthd <- "z"
    rval <- cbind(est, se, tval, pval)
    colnames(rval) <- cnames
    class(rval) <- "coeftest"
    attr(rval, "method") <- paste(mthd, "test of coefficients")
    return(rval)
}
## weights.stpm2 <- 
## function (object, ...) 
## {
##     wts <- object@weights
##     if (is.null(wts)) 
##         wts
##     else napredict(object@na.action, wts)
## }

require(sandwich)
coxph1 <- coxph(Surv(rectime,censrec==1)~hormon,data=brcancer)
update(coxph1,robust=TRUE)
sandwich(coxph1)
sandwich.stpm2(fit) # hurrah!


## require(lmtest)
## coeftest(coxph1)
## coeftest(coxph1,vcov.=sandwich(coxph1))
## coeftest(fit,sandwich(fit))


sandwich(fit)
sandwich(fit,bread.=bread.stpm2,meat.=meat.stpm2)


## some predictions
head(predict(fit,se.fit=TRUE,type="surv"))
head(predict(fit,se.fit=TRUE,type="hazard"))

## some plots
plot(fit,newdata=data.frame(hormon=0),type="hazard")

## time-varying coefficient
summary(fit.tvc <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,df=3,
                         tvc=list(hormon=3)))
anova(fit,fit.tvc) # compare with and without tvc

plot(fit.tvc,newdata=data.frame(hormon=0),type="hr",var="hormon")
                                        # no lines method: use add=TRUE
plot(fit.tvc,newdata=data.frame(hormon=1),type="hr",var="hormon",
     add=TRUE,ci=FALSE,line.col=2)

plot(fit.tvc,newdata=data.frame(hormon=0),type="sdiff",var="hormon")

plot(fit.tvc,newdata=data.frame(hormon=0),type="hdiff",var="hormon")

plot(fit.tvc,newdata=data.frame(hormon=0),type="hazard")
plot(fit.tvc,newdata=data.frame(hormon=1),type="hazard",line.col=2,ci=FALSE,add=TRUE)
## trace("predict", browser, exit=browser, signature = "stpm2")

set.seed(10101)
brcancer <- transform(brcancer, x=rlnorm(nrow(brcancer)))
summary(fit.tvc <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,df=3,
                     tvc.formula=~hormon:nsx(log(rectime),df=3)))



## cure model
## cf. http://www.pauldickman.com/survival/solutions/q37.do
### Data setup
require(foreign)
colon <- read.dta("http://www.pauldickman.com/survival/colon.dta")
popmort <- read.dta("http://www.pauldickman.com/survival/popmort.dta")
brcancer <- read.dta("http://www.stata-press.com/data/r11/brcancer.dta")
popmort <- transform(popmort, age=`_age`, year=`_year`, `_age`=NULL, `_year`=NULL)

save(colon,file="c:/usr/src/R/rstpm2/pkg/data/colon.rda")
save(popmort,file="c:/usr/src/R/rstpm2/pkg/data/popmort.rda")
save(brcancer,file="c:/usr/src/R/rstpm2/pkg/data/brcancer.rda")

## require(rstpm2)
popmort2 <- transform(popmort,exitage=age,exityear=year,age=NULL,year=NULL)
colon2 <- within(colon, {
  status <- ifelse(surv_mm>120.5,1,status)
  tm <- pmin(surv_mm,120.5)/12
  exit <- dx+tm*365.25
  sex <- as.numeric(sex)
  exitage <- pmin(floor(age+tm),99)
  exityear <- floor(yydx+tm)
})
colon2 <- merge(colon2,popmort2)

## compare relative survival without and with cure 
summary(fit0 <- stpm2(Surv(tm,status %in% 2:3)~I(year8594=="Diagnosed 85-94"),
                     data=colon2,
                     bhazard=colon2$rate, df=5)) ## CHECKED: same year8594 estimate as Stata
head(predict(fit0))
## estimate of failure at the end of follow-up
1-predict(fit0,data.frame(year8594 = unique(colon2$year8594),tm=max(colon2$tm)),type="surv",se.fit=TRUE)
plot(fit0,newdata=data.frame(year8594 = "Diagnosed 85-94"),ylim=0:1)
plot(fit0,newdata=data.frame(year8594 = "Diagnosed 75-84"),add=TRUE,line.col="red",rug=FALSE)
##
summary(fit <- stpm2(Surv(tm,status %in% 2:3)~I(year8594=="Diagnosed 85-94"),
                     data=colon2,
                     bhazard=colon2$rate,
                     df=5,cure=TRUE))
head(predict(fit))
## cure fractions (I need to add this to the predict function)
1-predict(fit,data.frame(year8594 = unique(colon2$year8594),tm=max(colon2$tm)),type="surv",se.fit=TRUE)
newdata1 <- data.frame(year8594 = "Diagnosed 85-94")
plot(fit,newdata=newdata1,add=TRUE,ci=FALSE,lty=2,rug=FALSE)
plot(fit,newdata=data.frame(year8594="Diagnosed 75-84"),add=TRUE,rug=FALSE,line.col="red",ci=FALSE,lty=2)

plot(fit,newdata=newdata1,type="hazard")
plot(fit,newdata=newdata1,type="cumhaz")


## http://www.pauldickman.com/survival/r/melanoma.relsurv.r
library(foreign)
library(survival)
library(relsurv)
# Download rates files from http://www.mortality.org/
# # 6. Life Tables By year of death (period) 1x1
# Save tables by gender in text files
# The transrate.hmd command translate these to R ratetables
Finlandpop <- transrate.hmd("c:/usr/tmp/mltper_1x1.txt","c:/usr/tmp/fltper_1x1.txt")

## The relsurv package requires time in days (exit and dx are dates of exit and diagnosis)
colon3 <- transform(colon2,tm.dd=as.numeric(exit-dx))
colon3$sex <- ifelse(colon2$sex==1,"male","female")
as.date <- function(x)
  if (inherits(x,"Date")) as.date(as.numeric(x)+3653) else date::as.date(x)
model1 <- rs.surv(Surv(tm.dd,status %in% 2:3)~year8594+ratetable(age=(X_age+0.5)*365.25,sex=sex,year=as.date(exit)),colon3,ratetable=Finlandpop)
plot(model1,lty=1:2)



                 
oldx <- 0:100
oldy <- (oldx-50)^2
oldy[c(20,30)] <- 0
old <- data.frame(x=oldx,y=oldy)
predict(lm(y~nsx(x,knots=c(25,50,75,95)),old)) # as per Stata
newx <- seq(min(oldx)/1.05,max(oldx)*1.05,length=101)
new <- data.frame(x=newx)
plot(oldx,oldy)
predict(lm(y~nsx(x,df=5,cure=TRUE),old))
sum(oldy)
terms(lm(y~nsx(x,df=5,cure=TRUE),old))
lm(y~nsx(x,df=5),old)


lines(newx,
      predict(lm(y~nsx(x,df=4,cure=FALSE),old),newdata=new),
      type="l") # oops
lines(newx,
      predict(lm(y~nsx(x,df=3),old),newdata=new),
      lty=2)


summary(fit <- stpm2(Surv(tm,status %in% 2:3)~I(year8594=="Diagnosed 85-94"),
                     data=colon2,
                     bhazard=colon2$rate,
                     logH.formula=~nsx(log(tm),df=6,stata=TRUE))) # okay
summary(fit <- stpm2(Surv(tm,status %in% 2:3)~I(year8594=="Diagnosed 85-94"),
                     data=colon2,
                     logH.formula=~nsx(log(tm),df=6,stata=TRUE))) # okay

## Stata
## stata.knots=c(4.276666164398193, 6.214608192443848, 6.7833251953125, 7.806289196014404)
stataKnots <- function(x,df) {
  intKnots <- round((1:(df-1))/df,2) # yes, Paul implicitly rounded to 2 dp
  logx <- log(x)
  c(min(logx),quantile(logx,intKnots,type=2),max(logx))
}
stata.knots <- stataKnots(subset(brcancer,censrec==1)$rectime,3)
## sapply(1:9,function(type) log(quantile(subset(brcancer,censrec==1)$rectime,c(0.33,0.67),type=type)))
summary(fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,
                     logH.args=list(knots=stata.knots[2:3],
                       Boundary.knots=stata.knots[c(1,4)])))
## formula specification for logH
summary(stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,
              logH.formula=~ns(log(rectime),df=3)))

pred <- predict(fit.tvc,newdata=data.frame(hormon=0:3),grid=TRUE,se.fit=TRUE,type="cumhaz")
pred.all <- cbind(pred,attr(pred,"newdata"))
## require(lattice)
## xyplot(Estimate ~ rectime, data=pred.all, group=hormon,type="l",xlab="Time")


## relative survival
brcancer <- transform(brcancer,rate0=10^(-5+x1/100))
summary(fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,bhazard=brcancer$rate0,df=3))
head(predict(fit,se.fit=TRUE))

## delayed entry
brcancer2 <- transform(brcancer,startTime=ifelse(hormon==0,rectime*0.5,0))
## debug(stpm2)
summary(fit <- stpm2(Surv(startTime,rectime,censrec==1)~hormon,data=brcancer2,
                     logH.formula=~nsx(log(rectime),df=3,stata=TRUE)))
head(predict(fit,se.fit=TRUE))
## delayed entry and tvc
summary(fit <- stpm2(Surv(startTime,rectime,censrec==1)~hormon,data=brcancer2,
                     tvc.formula=~hormon:nsx(log(rectime),df=3,stata=TRUE)))
head(predict(fit,se.fit=TRUE))



## multiple time scales
brcancer <- transform(brcancer,recyr=rectime/365.25)
## predictions from a simple model
summary(fit <- stpm2(Surv(recyr,censrec==1)~hormon+x1,data=brcancer,
                     logH.formula=~nsx(log(recyr),df=3,centre=log(50))))
head(predict(fit))
grid.x1 <- with(brcancer, seq(40,70,length=300))
newdata0 <- with(brcancer, data.frame(recyr=5,x1=grid.x1,hormon=0))
matplot(grid.x1,
        predict(fit,type="hr",newdata=newdata0,var="hormon",se.fit=TRUE), type="l")
## predictions with multiple time scales
summary(fit <- stpm2(Surv(recyr,censrec==1)~hormon,data=brcancer,
                     logH.formula=~nsx(log(recyr),df=3,centre=log(50)),
                     tvc.formula=~hormon:nsx(log(recyr+x1),df=2)))
matplot(grid.x1,
        predict(fit,type="hr",newdata=newdata0,var="hormon",se.fit=TRUE), type="l")





brcancer <- transform(brcancer,recyr=rectime/365.25,entry=recyr/2)
summary(fit <- stpm2(Surv(entry,recyr,censrec==1)~hormon,data=brcancer,
                     logH.formula=~nsx(log(recyr),df=3,centre=log(50)),
                     tvc.formula=~hormon:nsx(log(recyr+x1),df=2)))


summary(fit <- stpm2(Surv(recyr,censrec==1)~hormon+x1,data=brcancer,
                     logH.formula=~nsx(log(recyr),df=3,centre=log(50))))


plot(grid.x1,
     predict(fit,type="hr",newdata=newdata0,var="hormon",se.fit=TRUE)$fit, type="l")

plot(fit,newdata=data.frame(hormon=0,x1=50),var="hormon",type="hr")

head(predict(fit,type="hazard",newdata=newdata0))
head(predict(fit,type="hazard",newdata=transform(newdata0,hormon=1)))



newdata0 <- with(brcancer, data.frame(recyr=5+1,x1=grid.x1-1,hormon=0))
predict(fit,type="hr",newdata=newdata0,var="hormon")

summary(fit <- stpm2(Surv(recyr,censrec==1)~hormon+x1,data=brcancer,
                     logH.formula=~nsx(log(recyr),df=3,centre=log(50)),tvc=list(hormon=3)))

brcancer <- transform(brcancer, startAge=x1, endAge=x1+rectime/365)
summary(fit <- stpm2(Surv(startAge,endAge,censrec==1)~hormon,data=brcancer,
                     logH.formula=~nsx(log(endAge),df=3,centre=log(50)),tvc=list(hormon=3)))


## some simulated data: H_{weibull}(t)=(t/b)^a
n <- 1000
sim1 <- data.frame(age=seq(20,70,length=n),x=rep(0:1,each=n/2))
y <- rweibull(1000,shape=1,scale=1)



with(brcancer, plot(density(x1[censrec==1])))

summary(fit <- stpm2(Surv(recyr,censrec==1)~hormon,data=brcancer,logH.formula=~nsx(log(recyr),df=3,stata=TRUE)))




brcancer <- transform(brcancer,ageStart=rnorm(length(rectime),50,5))
brcancer <- transform(brcancer,ageStop=ageStart+rectime)
summary(fit <- stpm2(Surv(ageStart,ageStop,censrec==1)~hormon,data=brcancer,df=3))

brcancer3 <- transform(brcancer,startTime=ifelse(censrec==1,0,10))
summary(fit <- stpm2(Surv(startTime,rectime,censrec==1)~hormon,data=subset(brcancer,rectime>10),df=3))
summary(fit <- stpm2(Surv(startTime,rectime,censrec==1)~hormon,data=subset(brcancer3,rectime>10),df=3))

## check the performance time
refresh
require(rstpm2)
data(brcancer)
brcancer10 = do.call("rbind",lapply(1:10,function(i) brcancer))
system.time(summary(fit <- stpm2(Surv(rectime,censrec==1)~hormon,df=3,data=brcancer10)))
system.time(summary(fit <- stpm2(Surv(rectime,censrec==1)~hormon,df=3,data=brcancer10)))
system.time(summary(fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer10, logH.formula=~ns(log(rectime),df=4)+hormon:ns(log(rectime),df=3))))
system.time(summary(fit <- pstpm2(Surv(rectime,censrec==1)~hormon,data=brcancer10)))
system.time(summary(fit <- pstpm2(Surv(rectime,censrec==1)~1,data=brcancer10, smooth.formula=~s(log(rectime))+s(log(rectime),by=hormon))))

fit <- pstpm2(Surv(rectime,censrec==1)~hormon,data=brcancer10,trace=1)

fit <- pstpm2(Surv(rectime,censrec==1)~1,data=brcancer10, smooth.formula=~s(log(rectime))+s(log(rectime),by=hormon),trace=1,sp.init=c(1,1), reltol=list(outer=1e-5,search=1e-10,final=1e-10))

system.time(summary(fit <- pstpm2(Surv(rectime,censrec==1)~1,data=brcancer10, 
                                  smooth.formula=~s(log(rectime))+s(log(rectime),by=hormon),
                                  sp=c(0.006,0.0031),trace=1,outer_optim=2,criterion="GCV",
                                  reltol=list(outer=1e-5,search=1e-10,final=1e-10))))
## > fit@sp
## [1] 0.06104312 0.31430954

system.time(fit <- pstpm2(Surv(rectime,censrec==1)~1,data=brcancer10, smooth.formula=~s(log(rectime))+s(log(rectime),by=hormon),sp=c(1,1)))


nsx(1:10,df=3) - ns(1:10,df=3)
nsx(1:10,df=3,centre=3)
nsx(1:10,df=3,centre=3,Boundary.knots=c(2,8),derivs=c(1,1))
nsx(1:10,df=3,cure=TRUE)
nsxDeriv(1:10,df=3) - nsDeriv(1:10,df=3)
nsxDeriv(1:10,df=3,centre=5,derivs=c(1,1))
nsxDeriv(1:10,df=3,centre=5,cure=TRUE)

nsDeriv(1:10,df=3) - nsDeriv2(1:10,df=3)

## bug with calling mle2
require(bbmle)
mle2a <- function(...)
  mle2(...)
## some data
x <- 0:10
y <- c(26, 17, 13, 12, 20, 5, 9, 8, 5, 4, 8)
d <- data.frame(x,y)
## some fits
(fit0 <- mle2(y~dpois(lambda=ymean),start=list(ymean=mean(y)),data=d)) # okay
(fit0.2 <- mle2(y~dpois(lambda=ymean),start=list(ymean=mean(y)),data=d,
              control=list(parscale=2))) # okay
(fit1 <- mle2a(y~dpois(lambda=ymean),start=list(ymean=mean(y)),data=d)) # okay
(fit1.2 <- mle2a(y~dpois(lambda=ymean),start=list(ymean=mean(y)),data=d,
              control=list(parscale=2))) # FAILS



## stdReg::parfrailty documentation
library(stdReg)
library(survival)
     
## simulate data
n <- 1000
m <- 3
alpha <- 1.5
eta <- 1
phi <- 0.5
beta <- 1
id <- rep(1:n, each=m)
U <- rep(rgamma(n, shape=1/phi,scale=phi), each=m)
X <- rnorm(n*m)
## reparametrize scale as in rweibull function
weibull.scale <- alpha/(U*exp(beta*X))^(1/eta)
T <- rweibull(n*m, shape=eta, scale=weibull.scale)
## right censoring
C <- runif(n*m, 0,10)
D <- as.numeric(T<C)
T <- pmin(T, C)
## strong left-truncation
L <- runif(n*m, 0, 2)
incl <- T>L
incl <- ave(x=incl, id, FUN=sum)==m
dd <- data.frame(L, T, D, X, id)
dd <- dd[incl, ]  
##
fit <- parfrailty(formula=Surv(L, T, D)~X, data=dd, clusterid="id")
summary(fit)
##
library(rstpm2)
fit2 <- stpm2(formula=Surv(L, T, D)~X, data=dd, cluster=dd$id, smooth.formula=~log(T))
summary(fit2)

## ignore left truncation
fit <- parfrailty(formula=Surv(T, D)~X, data=dd, clusterid="id")
summary(fit)
fit2 <- stpm2(Surv(T, D)~X, data=dd, cluster=dd$id, smooth.formula=~log(T))
summary(fit2)
## normal random effect
fit2 <- stpm2(formula=Surv(T, D)~X, data=dd, cluster=dd$id, smooth.formula=~log(T), RandDist="LogN")
summary(fit2)



## end of examples ##




## ## * stata
## cd c:\Users\marcle\Documents\Home\
## clear
## webuse brcancer
## use brcancer
## stset rectime, f(censrec==1)
## cap program drop dopredictions
## program define dopredictions
##   preserve
##   predict hr, hrnumerator(hormon 1) ci
##   predict haz, hazard ci
##   predict surv, surv ci
##   predict sdiff, sdiff1(hormon 1) ci
##   list hr* in 1/5
##   list haz* surv* in 1/5
##   list sdiff* in 1/5
##   restore
## end

## * basic model
## stpm2 hormon, df(3) scale(h)
## dopredictions

## * cure 
## gen rate0=10^(-5+x1/100)
## stpm2 hormon, df(3) scale(h) cure bhazard(rate0)
## dopredictions

## * tvc
## stpm2 hormon, df(3) tvc(hormon) dftvc(3) scale(h)
## dopredictions

## * delayed entry
## preserve
##   replace _t0 = rectime*0.5 if hormon==0
##   stpm2 hormon, df(3) scale(h)
##   dopredictions
## restore

## * relative survival
## preserve  
##   gen rate0=10^(-5+x1/100)
##   stpm2 hormon, df(3) scale(h) bhazard(rate0)
##   dopredictions
## restore

## * test speed
## clear all
## set mem 100m
## use brcancer
## stset rectime, f(censrec==1)
## expand 100
## timer clear
## timer on 1
## stpm2 hormon, df(3) scale(h)
## timer off 1
## timer list
 

## hazard.pm = function(obj,tm,X,XD) # obj$par
## {
##   Xlocal=predict(X,newx=log(tm))
##   XDlocal=predict(XD,newx=log(tm))
##   with(obj,
##        c((XDlocal %*% par)/tm*exp(Xlocal %*% par)))
## }
## with(list(df=df,x=seq(0,3,length=100)[-1]),
##      {
##        plot(x,hazard.pm(fit,x,X,XD),type="l",ylim=c(0,2))
##        lines(x,dweibull(x,shape=1)/pweibull(x,shape=1,lower=FALSE),lty=2)
##      })
## ##
## require(deSolve)
## temp <- as.data.frame(ode(y=0,times=seq(0,10,length=100)[-1],
##                           func=function(t,state,parameters=NULL) list(exp(sin(2*pi*log(t))))))
## plot(temp,type="l")
## temp <- transform(temp, cum=`1`,logcum=log(`1`))
## with(temp,plot(log(time),logcum))
## temp1 <- temp[-1,]
## fit <- glm(log(cum)~log(time)+sin(2*pi*log(time))+cos(2*pi*log(time)),data=temp1)
## lines(log(temp1$time),predict(fit))
## ## In summary:
## ## we can model using sine and cosine terms for the log-cumulative hazard - for log(time).
