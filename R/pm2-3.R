## Utilities
## copied from stats:::format.perc
formating.perc <- 
    function (probs, digits) 
        paste(format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits), 
              "%")

## extension of ns() to include different boundary derivatives,
## centering and cure
nsx <- 
function (x, df = NULL, knots = NULL, intercept = FALSE,
          Boundary.knots = range(x),
          derivs = if (cure) c(2,1) else c(2,2),
          log=FALSE, # deprecated: only used in rstpm2:::stpm2Old
          centre = FALSE, cure = FALSE, stata.stpm2.compatible=FALSE) 
{
    nx <- names(x)
    x <- as.vector(x)
    nax <- is.na(x)
    if (nas <- any(nax)) 
        x <- x[!nax]
    if (!missing(Boundary.knots)) {
        Boundary.knots <- sort(Boundary.knots)
        outside <- (ol <- x < Boundary.knots[1L]) | (or <- x > 
            Boundary.knots[2L])
    }
    else outside <- FALSE
    if (!missing(df) && missing(knots)) {
        nIknots <- df - 1 - intercept + 4 - sum(derivs) 
        if (nIknots < 0) {
            nIknots <- 0
            warning("'df' was too small; have used ", 1 + intercept)
        }
        knots <- if (nIknots > 0) {
          knots <- if (!cure)
            seq.int(0, 1, length.out = nIknots + 2L)[-c(1L, 
                            nIknots + 2L)]
          else c(seq.int(0, 1, length.out = nIknots + 1L)[-c(1L, 
                                 nIknots + 1L)], 0.95)
          if (!stata.stpm2.compatible)
            stats::quantile(x[!outside], knots)
          else stats::quantile(x[!outside], round(knots,2), type=2)
        }
    }
    else nIknots <- length(knots)
    Aknots <- sort(c(rep(Boundary.knots, 4L), knots))
    if (any(outside)) {
        basis <- array(0, c(length(x), nIknots + 4L))
        if (any(ol)) {
            k.pivot <- Boundary.knots[1L]
            xl <- cbind(1, x[ol] - k.pivot)
            tt <- spline.des(Aknots, rep(k.pivot, 2L), 4, c(0, 
                1))$design
            basis[ol, ] <- xl %*% tt
        }
        if (any(or)) {
            k.pivot <- Boundary.knots[2L]
            xr <- cbind(1, x[or] - k.pivot)
            tt <- spline.des(Aknots, rep(k.pivot, 2L), 4, c(0, 
                1))$design
            basis[or, ] <- xr %*% tt
        }
        if (any(inside <- !outside)) 
            basis[inside, ] <- spline.des(Aknots, x[inside], 
                4)$design
    }
    else basis <- spline.des(Aknots, x, 4)$design
    const <- splineDesign(Aknots, rep(Boundary.knots, 3-derivs), 4, c(derivs[1]:2, derivs[2]:2))
    if (!intercept) {
        const <- const[, -1, drop = FALSE]
        basis <- basis[, -1, drop = FALSE]
    }
    qr.const <- qr(t(const))
    q.const <- qr.Q(qr.const, complete=TRUE)[, -(1L:(6-sum(derivs))), drop = FALSE] # NEW
    basis <- as.matrix((t(qr.qty(qr.const, t(basis))))[, -(1L:nrow(const)), drop = FALSE])
    n.col <- ncol(basis)
    if (nas) {
        nmat <- matrix(NA, length(nax), n.col)
        nmat[!nax, ] <- basis
        basis <- nmat
    }
    dimnames(basis) <- list(nx, 1L:n.col)
    if (is.numeric(centre)) {
      centreBasis <- nsx(centre,
                         knots=if (is.null(knots)) numeric(0) else knots,
                         Boundary.knots=Boundary.knots, 
                         intercept=intercept, derivs=derivs, centre=FALSE, log=log)
      oldAttributes <- attributes(basis)
      basis <- t(apply(basis,1,function(x) x-centreBasis))
      attributes(basis) <- oldAttributes
    }
    a <- list(degree = 3, knots = if (is.null(knots)) numeric(0) else knots, 
        Boundary.knots = Boundary.knots, intercept = intercept, derivs=derivs,
              centre=centre, log=log, q.const=q.const)
    attributes(basis) <- c(attributes(basis), a)
    class(basis) <- c("nsx", "basis", "matrix")
    basis
}
makepredictcall.nsx <- 
function (var, call) 
{
    if (as.character(call)[1L] != "nsx") 
        return(call)
    at <- attributes(var)[c("knots", "Boundary.knots", "intercept",
                            "derivs", "centre", "log")]
    xxx <- call[1L:2]
    xxx[names(at)] <- at
    xxx
}
predict.nsx <- 
function (object, newx, ...) 
{
    if (missing(newx)) 
        return(object)
    a <- c(list(x = newx), attributes(object)[c("knots", "Boundary.knots", 
        "intercept", "derivs", "centre", "log")])
    do.call("nsx", a)
}
Shat <- function(obj)
  {
    ## predicted survival for individuals (adjusted for covariates)
    newobj = survfit(obj,se.fit=FALSE)
    surv = newobj$surv
    rr = try(predict(obj,type="risk"),silent=TRUE)
    if (inherits(rr,"try-error"))
        rr <- 1
    surv2 = surv[match(obj$y[,ncol(obj$y)-1],newobj$time)]
    return(surv2^rr)
  }
replaceCall=function(obj,old,new) {
  if (is.atomic(obj) && length(obj)>1)
    return(as.call(c(quote(c),lapply(as.list(obj),replaceCall,old,new))))
  if (is.name(obj) || is.symbol(obj) || (is.atomic(obj) && length(obj)==1)) {
    if (obj==old) return(new)
    else return(obj)
  }
##   if (length(obj)==1 && length(obj[[1]])==1) {
##     if (obj==old) return(new)
##     else return(obj)
##   }
  as.call(lapply(obj,replaceCall,old,new))
}
replaceFormula=function(...) as.formula(replaceCall(...))
## replaceFormula(~f(a+b),quote(f),quote(g))
allCall=function(obj) {
  if (is.atomic(obj) && length(obj)==1) return(obj)
  if (is.atomic(obj) && length(obj)>1) return(as.call(c(quote(c),as.list(obj))))
  if (is.name(obj) || is.symbol(obj)) return(obj)
  as.call(lapply(obj,allCall))
}
## allCall(as.call(c(quote(ns),list(df=3,knots=c(1,2)))))[[2]]
vector2call=function(obj) {
  if (is.atomic(obj) && length(obj)==1) return(obj)
  if (is.atomic(obj) && length(obj)>1) return(as.call(c(quote(c),as.list(obj))))
  if (is.name(obj) || is.symbol(obj)) return(obj)
  lapply(obj,allCall) # is this correct?
}
## vector2call(list(df=3,knots=c(1,2)))
findSymbol <- function(obj,symbol) {
  if (is.symbol(obj) && obj==symbol) TRUE else
  if (is.symbol(obj)) FALSE else
  if (is.atomic(obj)) FALSE else
  Reduce(`|`,lapply(obj,findSymbol,symbol),FALSE)
}
rhs=function(formula) 
  if (length(formula)==3) formula[[3]] else formula[[2]]
lhs <- function(formula) 
  if (length(formula)==3) formula[[2]] else NULL
"rhs<-" = function(formula,value) {
  newformula <- formula
  newformula[[length(formula)]] <- value
  newformula
}
"lhs<-" <- function(formula,value) {
  if (length(formula)==2)
    as.formula(as.call(c(formula[[1]],value,formula[[2]])))
  else {
    newformula <- formula
    newformula[[2]] <- value
    newformula
  }
}

## numerically calculate the partial gradient \partial func_j \over \partial x_i
## (dim(grad(func,x)) == c(length(x),length(func(x)))
grad <- function(func,x,...) # would shadow numDeriv::grad()
  {
    h <- .Machine$double.eps^(1/3)*ifelse(abs(x)>1,abs(x),1)
    temp <- x+h
    h.hi <- temp-x
    temp <- x-h
    h.lo <- x-temp
    twoeps <- h.hi+h.lo
    nx <- length(x)
    ny <- length(func(x,...))
    if (ny==0L) stop("Length of function equals 0")
    df <- if(ny==1L) rep(NA, nx) else matrix(NA, nrow=nx,ncol=ny)
    for (i in 1L:nx) {
      hi <- lo <- x
      hi[i] <- x[i] + h.hi[i]
      lo[i] <- x[i] - h.lo[i]
      if (ny==1L)
        df[i] <- (func(hi, ...) - func(lo, ...))/twoeps[i]
      else df[i,] <- (func(hi, ...) - func(lo, ...))/twoeps[i]
      }
    return(df)
  }
## numerically calculate the gradient \partial func_i \over \partial x_i
## length(grad(func,x)) == length(func(x)) == length(x)
grad1 <- function(func,x,...,log.transform=FALSE)
{
    ny <- length(func(x,...))
    if (ny==0L) stop("Length of function equals 0")
    if (log.transform) {
        h <- .Machine$double.eps^(1/3)
        value <- (func(x*exp(h), ...) - func(x*exp(-h), ...))/2/h/x
    } else {
        h <- .Machine$double.eps^(1/3)*ifelse(abs(x)>1,abs(x),1)
        value <- (func(x+h, ...) - func(x-h, ...))/2/h
    }
    return(value)
}
## predict lpmatrix for an lm object
lpmatrix.lm <- 
  function (object, newdata, na.action = na.pass) {
    tt <- terms(object)
    if (!inherits(object, "lm")) 
      warning("calling predict.lm(<fake-lm-object>) ...")
    if (missing(newdata) || is.null(newdata)) {
      X <- model.matrix(object)
    }
    else {
      Terms <- delete.response(tt)
      m <- model.frame(Terms, newdata, na.action = na.action, 
                       xlev = object$xlevels)
      if (!is.null(cl <- attr(Terms, "dataClasses"))) 
        .checkMFClasses(cl, m)
      X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
    }
    X
  }
## fun: takes coef as its first argument
## requires: coef() and vcov() on the object
numDeltaMethod <- function(object, fun, gd=NULL,
                           conf.int=FALSE, level=0.95, ...) {
  coef <- coef(object)
  Sigma <- vcov(object)
  fit <- fun(coef,...)
  if (is.null(gd))
      gd <- grad(fun,coef,...)
  se.fit <- as.vector(sqrt(colSums(gd* (Sigma %*% gd))))
  if (!is.null(names(fit)))
      names(se.fit) <- names(fit)
  if(all(se.fit==0 | is.na(se.fit))) warning("Zero variance estimated. Do you need to pass a newdata argument to fun()?")
  df <- data.frame(fit = as.numeric(fit), se.fit = as.numeric(se.fit),
                   Estimate = as.numeric(fit), SE = as.numeric(se.fit))
  row.names(df) <- names(fit)
  if (conf.int) {
      a <- (1 - level)/2
      a <- cbind(a, 1 - a)
      fac <- qnorm(a)
      df$conf.low <-  df$fit+fac[,1]*df$se.fit
      df$conf.high <- df$fit+fac[,2]*df$se.fit
  }
  structure(df, # vcov=Sigma,
            class=c("predictnl","data.frame"))
}
"coef<-" <- function (x, value) 
    UseMethod("coef<-")
predictnl <- function (object, ...) 
  UseMethod("predictnl")
setGeneric("predictnl")
"coef<-.default" <- function(x,value) {
    x$coefficients <- value
    x
}
"coef<-.mle2" <- function(x,value) {
    x@fullcoef <- value
    x
}
predictnl.default <- function(object,fun,newdata=NULL,gd=NULL,...)
  {
      if (!is.null(newdata) || "newdata" %in% names(formals(fun))) {
          local1 <- function(coef,newdata,...)
              {
                  coef(object) <- coef
                  fun(object,newdata=newdata,...)
              }
          numDeltaMethod(object, local1, newdata=newdata, gd=gd, ...)
      }
      else {
          local2 <- function(coef,...)
              {
                  coef(object) <- coef
                  fun(object,...)
              }
          numDeltaMethod(object, local2, gd=gd, ...)
      }
  }
setMethod("predictnl", "mle2", function(object,fun,newdata=NULL,gd=NULL,...)
  {
    if (is.null(newdata) && !is.null(object@data))
      newdata <- object@data
    localf <- function(coef,...)
      {
        object@fullcoef <- coef # changed from predictnl.default()
        fun(object,...)
      }
    numDeltaMethod(object,localf,newdata=newdata,gd=gd,...)
  })
confint.predictnl <- function(object,parm,level=0.95,...) {
    cf <- object$fit
    pnames <- names(cf)
    if (is.null(pnames))
        pnames <- 1:length(cf)
    if (missing(parm)) 
        parm <- pnames
    else if (is.numeric(parm)) 
        parm <- pnames[parm]
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    pct <- formating.perc(a, 3)
    fac <- qnorm(a)
    ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, 
        pct))
    ses <- object$se.fit[parm]
    ci[] <- as.vector(cf[parm]) + ses %o% fac
    ci
}
predictnl.lm <- 
function (object, fun, newdata = NULL, ...) 
{
    if (is.null(newdata) && "newdata" %in% names(formals(fun))) {
        stopifnot(!is.null(object$data))
        newdata <- object$data
    }
    predictnl.default(object, fun, newdata, ...)
}
## setMethod("predictnl", "mle", function(object,fun,gd=NULL,...)
##   {
##     localf <- function(coef,...)
##       {
##         object@fullcoef = coef # changed from predictnl.default()
##         fun(object,...)
##       }
##     numDeltaMethod(object,localf,gd=gd,...)
##   })
predict.formula <- function(object,data,newdata,na.action,type="model.matrix",
                            ...) 
{
  mf <- match.call(expand.dots = FALSE)
  type <- match.arg(type)
  m <- match(c("formula", "data", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  xlevels <-.getXlevels(mt, mf)
  mfnew <- model.frame(mt, newdata, na.action=na.action, xlev=xlevels)
  if (!is.null(cl <- attr(mt, "dataClasses"))) .checkMFClasses(cl, mfnew)
  model.matrix(mt, mfnew, contrasts=contrasts)
}
`%call+%` <- function(left,right) call("+",left,right)
`%call-%` <- function(left,right) call("-",left,right)
##
bread.stpm2 <- function (x, ...) {
  rval <- vcov(x) * nrow(x@y)
  dimnames(rval) <- list(names(coef(x)), names(coef(x)))
  return(rval)
}
estfun.stpm2 <- function(obj, weighted=FALSE, ...) {
  rr <- t(grad(obj@logli,coef(obj)))
  colnames(rr) <- names(coef(obj))
  if (weighted)
    rr <- rr * obj@weights
  rr
}
applyTapplySum <- function(X,index) apply(X, 2, function(col) tapply(col, index, sum))
meat.stpm2 <- function (x, adjust = FALSE, cluster=NULL, ...) 
{
    psi <- estfun.stpm2(x, ...)
    k <- NCOL(psi)
    n <- NROW(psi)
    if (!is.null(cluster))
        psi <- applyTapplySum(as.matrix(psi),cluster)
    rval <- crossprod(as.matrix(psi))/n
    if (adjust) 
        rval <- n/(n - k) * rval
    rownames(rval) <- colnames(rval) <- colnames(psi)
    return(rval)
}
sandwich.stpm2 <- 
function (x, bread. = bread.stpm2, meat. = meat.stpm2, ...) 
{
    if (is.function(bread.)) 
        bread. <- bread.(x)
    if (is.function(meat.)) 
        meat. <- meat.(x, ...)
    n <- NROW(estfun.stpm2(x))
    return(1/n * (bread. %*% meat. %*% bread.))
}
incrVar <- function(var,increment=1) {
  ##var <- deparse(substitute(var))
  ##function(data) "$<-"(data,var,"$"(data,var)+increment) # FAILS
  n <- length(var)
  if (n>1 && length(increment)==1)
    increment <- rep(increment,n)
  function(data) {
    for (i in 1:n) {
      data[[var[i]]] <- data[[var[i]]] + increment[i]
    }
    data
  }
}
cloglog <- function(x) log(-log(x))
cexpexp <- function(x) exp(-exp(x))
setOldClass("terms")
setClassUnion("listOrNULL",c("list","NULL"))
setClassUnion("nameOrcall",c("name","call"))
setClassUnion("nameOrcallOrNULL",c("name","call","NULL"))
##setClassUnion("numericOrNULL",c("numeric","NULL"))
setOldClass("Surv")
setOldClass("lm")
expit <- function(x) {
    ifelse(x==-Inf, 0, ifelse(x==Inf, 1, 1/(1+exp(-x))))
}
logit <- function(p) {
    ifelse(p==0, -Inf, ifelse(p==1, Inf, log(p/(1-p))))
} # numerical safety for large values?
## check: weights
##
## adapted from ordinal::drop.coef
which.dim <- function (X, silent = TRUE) 
{
    stopifnot(is.matrix(X))
    silent <- as.logical(silent)[1]
    qr.X <- qr(X, tol = 1e-07, LAPACK = FALSE)
    if (qr.X$rank == ncol(X)) 
        return(TRUE)
    if (!silent) 
        message(gettextf("design is column rank deficient so dropping %d coef", 
            ncol(X) - qr.X$rank))
    return(qr.X$pivot[1:qr.X$rank])
}

## link families
link.PH <- list(link=function(S) log(-log(as.vector(S))),
                ilink=function(eta) exp(-exp(as.vector(eta))),
                gradS=function(eta,X) -exp(as.vector(eta))*exp(-exp(as.vector(eta)))*X,
                h=function(eta,etaD) as.vector(etaD)*exp(as.vector(eta)),
                H=function(eta) exp(as.vector(eta)),
                gradh=function(eta,etaD,obj) obj$XD*exp(as.vector(eta))+obj$X*as.vector(etaD)*exp(as.vector(eta)),
                gradH=function(eta,obj) obj$X*exp(as.vector(eta)))
link.PO <- list(link=function(S) -logit(as.vector(S)),
                ilink=function(eta) expit(-as.vector(eta)),
                gradS=function(eta,X) -(exp(as.vector(eta))/(1+exp(as.vector(eta)))^2)*X,
                H=function(eta) log(1+exp(as.vector(eta))),
                h=function(eta,etaD) as.vector(etaD)*exp(as.vector(eta))*expit(-as.vector(eta)),
                gradh=function(eta,etaD,obj) {
                    as.vector(etaD)*exp(as.vector(eta))*obj$X*expit(-as.vector(eta)) -
                        exp(2*as.vector(eta))*obj$X*as.vector(etaD)*expit(-as.vector(eta))^2 +
                            exp(as.vector(eta))*obj$XD*expit(-as.vector(eta))
                    },
                gradH=function(eta,obj) obj$X*exp(as.vector(eta))*expit(-as.vector(eta)))
link.probit <-
    list(link=function(S) -qnorm(as.vector(S)),
         ilink=function(eta) pnorm(-as.vector(eta)),
         gradS=function(eta,X) -dnorm(-as.vector(eta))*X,
         H=function(eta) -log(pnorm(-as.vector(eta))),
         h=function(eta,etaD) dnorm(as.vector(eta))/pnorm(-as.vector(eta))*as.vector(etaD),
         gradh=function(eta,etaD,obj) {
             -as.vector(eta)*obj$X*dnorm(as.vector(eta))*as.vector(etaD)/pnorm(-as.vector(eta)) +
                 obj$X*dnorm(as.vector(eta))^2/pnorm(-as.vector(eta))^2*as.vector(etaD) +
                     dnorm(as.vector(eta))/pnorm(-as.vector(eta))*obj$XD
         },
         gradH=function(eta,obj) obj$X*dnorm(as.vector(eta))/pnorm(-as.vector(eta)))
link.AH <- list(link=function(S) -log(S),
                ilink=function(eta) exp(-as.vector(eta)),
                gradS=function(eta,X) -as.vector(exp(-as.vector(eta)))*X,
                h=function(eta,etaD) as.vector(etaD),
                H=function(eta) as.vector(eta),
                gradh=function(eta,etaD,obj) obj$XD,
                gradH=function(eta,obj) obj$X)
link.AO <- function(theta) { # Aranda-Ordaz
    if (theta==0) {
        return(link.PH)
        } else {
            list(link = function(S) log((S^(-theta)-1)/theta),
                 ilink = function(eta) exp(-log(theta*exp(as.vector(eta))+1)/theta),
                 gradS = function(eta,X) -as.vector(exp(as.vector(eta))*exp(-log(theta*exp(as.vector(eta))+1)/theta)/(1+theta*exp(as.vector(eta))))*X,
                 H = function(eta) log(theta*exp(as.vector(eta))+1)/theta,
                 h = function(eta,etaD) exp(as.vector(eta))*as.vector(etaD)/(theta*exp(as.vector(eta))+1),
                 gradH = function(eta,obj) exp(as.vector(eta))*obj$X/(1+theta*exp(as.vector(eta))),
                 gradh = function(eta,etaD,obj) {
                     eta <- as.vector(eta)
                     etaD <- as.vector(etaD)
                     ((theta*exp(2*eta)+exp(eta))*obj$XD+exp(eta)*etaD*obj$X) /
                         (theta*exp(eta)+1)^2
                 })
        }
    }
## fd <- function(f,x,eps=1e-5) (f(x+eps)-f(x-eps))/2/eps
fd <- function(f,x,eps=1e-5)
    t(sapply(1:length(x),
             function(i) {
                 upper <- lower <- x
                 upper[i]=x[i]+eps
                 lower[i]=x[i]-eps
                 (f(upper)-f(lower))/2/eps
             }))
## test code for the link functions
if (FALSE) {
    Xstar <- cbind(1,1:3) # beta[1] + beta[2]*t
    betastar <- c(-4,0.5)
    XDstar <- cbind(0,Xstar[,2])
    etastar <- as.vector(Xstar %*% betastar)
    etaDstar <- as.vector(XDstar %*% betastar)
    obj <- list(X=Xstar,XD=XDstar)
    for(link in list(rstpm2:::link.PH,rstpm2:::link.PO,rstpm2:::link.probit,rstpm2:::link.AH,rstpm2:::link.AO(.3))) {
        print(rstpm2:::fd(function(beta) link$ilink(Xstar%*%beta), betastar)-t(link$gradS(etastar,Xstar)))
        print(rstpm2:::fd(function(beta) link$h(Xstar%*%beta, XDstar%*%beta), betastar)-t(link$gradh(etastar,etaDstar,obj)))
        print(rstpm2:::fd(function(beta) link$H(Xstar%*%beta), betastar)-t(link$gradH(etastar,obj)))
    }
}
bhazard <- function(x) x

## general link functions
setClass("stpm2", representation(xlevels="list",
                                 contrasts="listOrNULL",
                                 terms="terms",
                                 logli="function",
                                 ## weights="numericOrNULL",
                                 lm="lm",
                                 timeVar="character",
                                 time0Var="character",
                                 timeExpr="nameOrcall",
                                 time0Expr="nameOrcallOrNULL",
                                 delayed="logical",
                                 interval="logical",
                                 frailty="logical",
                                 model.frame="list",
                                 call.formula="formula",
                                 x="matrix",
                                 xd="matrix",
                                 termsd="terms",
                                 Call="call",
                                 y="Surv",
                                 link="list",
                                 args="list"
                                 ),
         contains="mle2")

gsm.control <- function(parscale=1,
                               maxit=300,
                               optimiser=c("BFGS","NelderMead"),
                               trace=0,
                               nodes=9,
                               adaptive=TRUE,
                               kappa.init=1,
                               maxkappa=1e3,
                               suppressWarnings.coxph.frailty=TRUE,
                               robust_initial=FALSE,
                        bhazinit=0.1,
                        eps.init=1e-5,
                               use.gr=TRUE,
                               penalty=c("logH","h"),
                               outer_optim=1,
                               reltol.search=1e-10, reltol.final=1e-10, reltol.outer=1e-5,
                               criterion=c("GCV","BIC")) {
    stopifnot.logical <- function(arg)
        stopifnot(is.logical(arg) || (is.numeric(arg) && arg>=0))
    stopifnot(parscale>0)
    stopifnot(maxit>1)
    optimiser <- match.arg(optimiser)
    stopifnot(trace>=0)
    stopifnot(nodes>=3)
    stopifnot.logical(adaptive)
    stopifnot(maxkappa>0)
    stopifnot.logical(suppressWarnings.coxph.frailty)
    stopifnot.logical(robust_initial)
    stopifnot(bhazinit>0)
    stopifnot.logical(use.gr)
    stopifnot(kappa.init>0)
    penalty <- match.arg(penalty)
    stopifnot(outer_optim %in% c(0,1))
    stopifnot(reltol.search>0)
    stopifnot(reltol.final>0)
    stopifnot(reltol.outer>0)
    criterion <- match.arg(criterion)
    list(mle2.control=list(parscale=parscale, maxit=maxit),
         optimiser=optimiser, trace=trace, nodes=nodes,
         adaptive=adaptive, maxkappa=maxkappa,
         suppressWarnings.coxph.frailty=suppressWarnings.coxph.frailty,
         robust_initial=robust_initial, bhazinit=bhazinit, eps.init=eps.init,
         use.gr=use.gr,
         kappa.init=kappa.init, penalty=penalty, outer_optim=outer_optim,
         reltol.search=reltol.search, reltol.final=reltol.final,
         reltol.outer=reltol.outer,
         criterion=criterion)
}

gsm <- function(formula, data, smooth.formula = NULL, smooth.args = NULL,
                df = 3, cure = FALSE,
                tvc = NULL, tvc.formula = NULL,
                control = list(), init = NULL,
                weights = NULL, robust = FALSE, baseoff = FALSE,
                timeVar = "", time0Var = "", use.gr = NULL,
                optimiser=NULL, log.time.transform=TRUE,
                reltol=NULL, trace = NULL,
                link.type=c("PH","PO","probit","AH","AO"), theta.AO=0,
                contrasts = NULL, subset = NULL,
                robust_initial=NULL,
                ## arguments for specifying the Cox model for initial values (seldom used)
                coxph.strata = NULL, coxph.formula = NULL,
                ## deprecated arguments (for removal)
                logH.formula = NULL, logH.args = NULL,
                ## arguments specific to relative survival
                bhazard = NULL, bhazinit=NULL,
                ## arguments specific to the fraility models and copulas
                copula=FALSE,
                frailty = !is.null(cluster) & !robust & !copula, cluster = NULL, logtheta=NULL,
                nodes=NULL, RandDist=c("Gamma","LogN"), recurrent = FALSE,
                adaptive = NULL, maxkappa = NULL,
                ## arguments specific to the penalised models
                sp=NULL, criterion=NULL, penalty=NULL,
                smoother.parameters=NULL, Z=~1, outer_optim=NULL,
                alpha=1, sp.init=1,
                penalised=FALSE,
                ## other arguments
                ...) {
    link.type <- match.arg(link.type)
    link <- switch(link.type,PH=link.PH,PO=link.PO,probit=link.probit,AH=link.AH,
                   AO=link.AO(theta.AO))
    RandDist <- match.arg(RandDist)
    ## handle deprecated args
    deprecated.arg <- function(name)
        if (!is.null(val <- get(name))) {
            warning(sprintf("%s argument is deprecated. Use control=list(%s=value)",
                            name,name))
            control[[name]] <- val
            TRUE
        } else FALSE
    deprecated <- lapply(c("optimiser","reltol","trace","nodes","adaptive","maxkappa",
                           "robust_initial","bhazinit", "use.gr", "penalty", "outer_optim",
                           "criterion"),
                         deprecated.arg)
    if(!is.null(reltol)) {
        warning("reltol argument is deprecated.\nUse control=list(reltol.search=value,reltol.final=value,reltol.outer=value).")
        stopifnot(all(c("search","final","outer") %in% names(reltol)))
        control$reltol.search <- reltol$search
        control$reltol.final <- reltol$final
        control$reltol.outer <- reltol$outer
    }
    ## read from control list
    control <- do.call(gsm.control, control)
    ## logH.formula and logH.args are DEPRECATED
    if (is.null(smooth.formula) && !is.null(logH.formula)) {
        warning("logH.formula is deprecated - use smooth.formula")
        smooth.formula <- logH.formula
    }
    if (is.null(smooth.args) && !is.null(logH.args)) {
        warning("logH.args is deprecated - use smooth.formula")
        smooth.args <- logH.args
    }
    ## set na.action=na.pass (and reset at end)
    na.action.old <- options()[["na.action"]]
    options(na.action = "na.pass")
    on.exit(options(na.action = na.action.old))
    if (copula && control$use.gr) {
        control$use.gr <- FALSE
    }
    ##
    ## set up the data
    ## ensure that data is a data frame
    temp.formula <- formula
    if (!is.null(smooth.formula)) rhs(temp.formula) <-rhs(temp.formula) %call+% rhs(smooth.formula)
    raw.data <- data
    ## data <- get_all_vars(temp.formula, raw.data)
    ## parse the function call
    Call <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "contrasts", "weights"),
               names(mf), 0L)
    mf <- mf[c(1L, m)]
    ##
    ## parse the event expression
    eventInstance <- eval(lhs(formula),envir=data)
    stopifnot(length(lhs(formula))>=2)
    delayed <- length(lhs(formula))>=4 # indicator for multiple times (cf. strictly delayed)
    surv.type <- attr(eventInstance,"type")
    if (surv.type %in% c("interval2","left","mstate"))
        stop("gsm not implemented for Surv type ",surv.type,".")
    counting <- attr(eventInstance,"type") == "counting"
    interval <- attr(eventInstance,"type") == "interval"
    timeExpr <- lhs(formula)[[if (delayed && !interval) 3 else 2]] # expression
    eventExpr <- if (interval) lhs(formula)[[4]] else lhs(formula)[[if (delayed) 4 else 3]]
    if (interval)
        time2Expr <- lhs(formula)[[3]]
    if (timeVar == "")
        timeVar <- all.vars(timeExpr)
    ##
    ## set up the formulae
    if (is.null(smooth.formula) && is.null(smooth.args)) {
        if (penalised) {
            smooth.args$k <- -1
        }
        else {
            smooth.args$df <- df
            if (cure) smooth.args$cure <- cure
        }
    }
    smoother <- if (penalised) quote(s) else quote(nsx)
    if (is.null(smooth.formula)) {
        smooth.formula <- as.formula(call("~",as.call(c(smoother,call("log",timeExpr),
                                                        vector2call(smooth.args)))))
        if (link.type=="AH") {
            if (penalised) {
                interaction <- function(expr) {
                    if (is.name(expr))
                        call(":",expr,timeExpr)
                    else if(is.name(expr[[1]]) && as.character(expr[[1]])=="+")
                        interaction(expr[[2]]) %call+% interaction(expr[[3]])
                    else call(":",expr,timeExpr)
                }
                ## interaction(rstpm2:::rhs(~a+I(b)+c-1)) # error
                smooth.formula <-
                    as.formula(call("~",
                                    interaction(rhs(formula)) %call+%
                                    as.call(c(smoother,timeExpr, vector2call(smooth.args)))))
            }
            else {
                smooth.formula <-
                    as.formula(call("~",as.call(c(smoother,timeExpr,
                                                  vector2call(smooth.args))) %call+%
                                        as.call(c(as.name(":"),rhs(formula),timeExpr))))
            }
        }
    }
    if (!penalised && is.null(tvc.formula) && !is.null(tvc)) {
        tvc.formulas <-
            lapply(names(tvc), function(name)
                call(":",
                     call("as.numeric",as.name(name)), # is this a good idea?
                     as.call(c(smoother,
                               if (link.type=="AH") timeExpr else call("log",timeExpr),
                               if (penalised) vector2call(list(k=tvc[[name]]))
                               else vector2call(if (cure) list(cure=cure,df=tvc[[name]])
                                                else list(df=tvc[[name]]))))))
        if (length(tvc.formulas)>1)
            tvc.formulas <- list(Reduce(`%call+%`, tvc.formulas))
        tvc.formula <- as.formula(call("~",tvc.formulas[[1]]))
    }
    if (penalised && is.null(tvc.formula) && !is.null(tvc)) {
        tvc.formulas <-
            lapply(names(tvc), function(name)
                as.call(c(quote(s),
                          call("log",timeExpr),
                          vector2call(list(by=as.name(name),k=tvc[[name]])))))
        if (length(tvc.formulas)>1)
            tvc.formulas <- list(Reduce(`%call+%`, tvc.formulas))
        tvc.formula <- as.formula(call("~",tvc.formulas[[1]]))
        ## remove any main effects (special case: a single term)
        Terms <- terms(formula)
        constantOnly <- FALSE
        for (name in names(tvc)) {
            .i <- match(name,attr(Terms,"term.labels"))
            if (!is.na(.i)) {
                if (.i==1 && length(attr(Terms,"term.labels"))==1) {
                    constantOnly <- TRUE
                } else Terms <- stats::drop.terms(Terms,.i,keep=TRUE)
            }
        }
        formula <- formula(Terms)
        if (constantOnly)
            rhs(formula) <- quote(1)
        rm(constantOnly,Terms)
    }
    if (!is.null(tvc.formula)) {
        rhs(smooth.formula) <- rhs(smooth.formula) %call+% rhs(tvc.formula)
    }
    if (baseoff)
        rhs(smooth.formula) <- rhs(tvc.formula)
    full.formula <- formula
    if (link.type == "AH") {
        rhs(full.formula) <- rhs(smooth.formula)
    } else {
        rhs(full.formula) <- rhs(formula) %call+% rhs(smooth.formula)
    }
    Z.formula <- Z
    ##
    tf <- terms.formula(smooth.formula, specials = c("s", "te"))
    terms <- attr(tf, "term.labels")
    ##
    ##
    ## Different models:
    ## lm (full.formula - cluster), survreg (formula - cluster), coxph (formula - cluster), full
    lm.formula <- formula(full.formula)
    base.formula <- formula
    ## Specials:
    specials.names <- c("cluster","bhazard","offset")
    specials <- attr(terms.formula(formula, specials.names), "specials")
    spcall <- mf
    spcall[[1]] <- quote(stats::model.frame)
    spcall$formula <- terms(formula, specials.names, data = data)
    mf2 <- eval(spcall, parent.frame())
    if (any(!sapply(specials,is.null))) {
        cluster.index <- specials$cluster
        bhazard.index <- specials$bhazard
        if (length(cluster.index)>0) {
            cluster <- mf2[, cluster.index]
            frailty = !is.null(cluster) && !robust && !copula
            cluster.index2 <- attr(terms.formula(full.formula, "cluster"), "specials")$cluster
            base.formula <- formula(stats::drop.terms(terms(mf2), cluster.index - 1, keep.response=TRUE))
            lm.formula <- formula(stats::drop.terms(terms(full.formula), cluster.index2 - 1))
        } else {
            cluster.index2 = NULL
        }
        if (length(bhazard.index)>0) {
            bhazard <- mf2[, bhazard.index]
            bhazard.index2 <- attr(terms.formula(full.formula, "bhazard"), "specials")$bhazard
            termobj = terms(mf2)
            dropped.terms = if(length(attr(termobj,"term.labels"))==1) reformulate("1", response=termobj[[2L]], intercept=TRUE, env=environment(termobj)) else stats::drop.terms(terms(mf2), bhazard.index - 1, keep.response=TRUE)
            base.formula <- formula(dropped.terms)
            lm.formula <- formula(stats::drop.terms(terms(full.formula), bhazard.index2 - 1))
        } else {
            bhazard.index2 = NULL
        }
        if (length(cluster.index)>0 && length(bhazard.index)>0) {
            dropped.terms = if(length(attr(termobj,"term.labels"))==2) reformulate("1", response=termobj[[2L]], intercept=TRUE, env=environment(termobj)) else stats::drop.terms(terms(mf2), c(cluster.index,bhazard.index) - 1, keep.response=TRUE)
            base.formula <- formula(dropped.terms)
            ## base.formula <- formula(stats::drop.terms(terms(mf2),
            ##                                           c(cluster.index,bhazard.index) - 1,
            ##                                           keep.response=TRUE))
            lm.formula <- formula(stats::drop.terms(terms(full.formula),
                                                    c(cluster.index2,bhazard.index2) - 1))
        }
        ## rm(mf2,spcall)
    }
    ## offset
    offset <- as.vector(stats::model.offset(mf2))
    rm(mf2)
    ## deprecated code for cluster
    cluster <- substitute(cluster)
    cluster <- switch(class(cluster),
                      integer=cluster,
                      numeric=cluster,
                      call=eval(cluster, data, parent.frame()),
                      NULL=NULL,
                      name=eval(cluster, data, parent.frame()),
                      character=data[[cluster]])
    ## deprecated code for bhazard
    bhazard <- substitute(bhazard)
    bhazard <- switch(class(bhazard),
                      integer=bhazard,
                      numeric=bhazard,
                      call=eval(bhazard,data,parent.frame()),
                      NULL=rep(0,nrow(data)),
                      name=eval(bhazard, data, parent.frame()),
                      character=data[[bhazard]])
    ##
    subset.expr <- substitute(subset)
    if(inherits(subset.expr,"NULL")) subset.expr <- TRUE
    .include <- complete.cases(model.matrix(formula, data)) &
        !is.na(eval(eventExpr,data,parent.frame())) &
        eval(subset.expr,data,parent.frame())
    options(na.action = na.action.old)
    if (!interval) {
        time <- eval(timeExpr,data,parent.frame())
        if (any(is.na(time))) warning("Some event times are NA")
        if (any(ifelse(is.na(time),FALSE,time<=0))) warning("Some event times <= 0")
        .include <- .include & ifelse(is.na(time), FALSE, time>0)
    }
    time0Expr <- NULL # initialise
    if (delayed) {
        time0Expr <- lhs(formula)[[2]]
        if (time0Var == "")
            time0Var <- all.vars(time0Expr)
        time0 <- eval(time0Expr,data,parent.frame())
        if (any(is.na(time0))) warning("Some entry times are NA")
        if (any(ifelse(is.na(time0),FALSE,time0<0))) warning("Some entry times < 0")
        .include <- .include & ifelse(is.na(time0), FALSE, time0>=0)
    }
    if (!is.null(substitute(weights)))
        .include <- .include & !is.na(eval(substitute(weights),data,parent.frame()))
    ##
    if (!is.null(cluster))
        .include <- .include & !is.na(cluster)
    .include <- .include & !is.na(bhazard)
    if (!is.null(offset))
        .include <- .include & !is.na(offset)
    excess <- !all(bhazard==0)
    ##
    data <- data[.include, , drop=FALSE]
    ## we can now evaluate over data
    time <- eval(timeExpr, data, parent.frame())
    if (delayed) {
        time0 <- eval(time0Expr, data, parent.frame())
    }
    event <- eval(eventExpr,data,parent.frame())
    if (!interval)
        event <- if (length(unique(event))==1) rep(TRUE, length(event))
                 else event <- event > min(event)
    nevent <- sum(event)
    if (!is.null(cluster))
        cluster <- cluster[.include]
    bhazard <- bhazard[.include]
    offset <- if (is.null(offset)) rep(0,sum(.include)) else offset[.include]
    if (copula && delayed)
        stop("Copulas not currently defined for delayed entry")
    ## setup for initial values
    if (!interval) {
        ## Cox regression
        coxph.call <- mf
        coxph.call[[1L]] <- quote(survival::coxph)
        coxph.strata <- substitute(coxph.strata)
        coxph.call$data <- quote(coxph.data)
        coxph.data <- data # ?
        coxph.call$formula <- base.formula
        if (!is.null(coxph.formula)) {
            temp <- coxph.call$formula
            rhs(temp) <- rhs(temp) %call+% rhs(coxph.formula)
            coxph.call$formula <- temp
        }
        if (!is.null(coxph.strata)) {
            temp <- coxph.call$formula
            rhs(temp) <- rhs(temp) %call+% call("strata",coxph.strata)
            coxph.call$formula <- temp
        }
        coxph.call$model <- TRUE
        coxph.obj <- eval(coxph.call, coxph.data)
        y <- model.extract(model.frame(coxph.obj),"response")
        data$logHhat <- if (is.null(bhazard)) {
                            link$link(pmin(1-control$eps.init,
                                           pmax(control$eps.init,Shat(coxph.obj))))
                        } else  link$link(pmin(1-control$eps.init,
                                               pmax(control$eps.init,Shat(coxph.obj)/
                                                                     exp(-control$bhazinit*bhazard*time))))
        if ((frailty || copula) && is.null(logtheta)) { # fudge for copulas
            coxph.data$.cluster <- as.vector(unclass(factor(cluster)))
            coxph.formula <- coxph.call$formula
            rhs(coxph.formula) <- rhs(coxph.formula) %call+%
                call("frailty",as.name(".cluster"),
                     distribution=switch(RandDist,LogN="gaussian",Gamma="gamma"))
            coxph.call$formula <- coxph.formula
            coxph.obj <- if (control$suppressWarnings.coxph.frailty) suppressWarnings(eval(coxph.call, coxph.data)) else eval(coxph.call, coxph.data)
            logtheta <- log(coxph.obj$history[[1]]$theta)
        }
    }
    if (interval) {
        ## survref regression
        survreg.call <- mf
        survreg.call[[1L]] <- as.name("survreg")
        survreg.call$formula <- base.formula
        survreg.obj <- eval(survreg.call, envir=parent.frame())
        weibullShape <- 1/survreg.obj$scale
        weibullScale <- predict(survreg.obj)
        y <- model.extract(model.frame(survreg.obj),"response")
        data$logHhat <- link$link(pmin(1-control$eps.init,
                                       pmax(control$eps.init,
                                            pweibull(time,weibullShape,weibullScale,lower.tail=FALSE))))
        ##
        if ((frailty || copula) && is.null(logtheta)) {
            logtheta <- -1
        }
    }
    ## initial values and object for lpmatrix predictions
    lm.call <- mf
    lm.call[[1L]] <- if (penalised) quote(gam) else quote(stats::lm)
    lhs(lm.formula) <- quote(logHhat) # new response
    lm.call$formula <- lm.formula
    dataEvents <- data[event,]
    if (interval) {
        dataEvents <- data[event>0, , drop=FALSE]
    }
    lm.call$data <- quote(dataEvents) # events only
    if (penalised) {
        lm.call$sp <- sp
        if (is.null(sp) && !is.null(sp.init) && (length(sp.init)>1 || sp.init!=1))
            lm.call$sp <- sp.init
    }
    lm.obj <- eval(lm.call)
    ## re-run gam if sp.init==1 (default)
    if (penalised && is.null(sp) && !is.null(sp.init) && length(sp.init)==1 && sp.init==1) {
        sp.init <- lm.call$sp <- rep(sp.init,length=length(lm.obj$sp))
        lm.obj <- eval(lm.call)
    }
    has.init <- !is.null(init)
    if (!has.init) {
        init <- coef(lm.obj)
    } else {
        stopifnot(length(init) == length(coef(lm.obj)))
    }
    ##
    ## set up mf and wt
    mt <- terms(lm.obj)
    mf <- model.frame(lm.obj)
    ## wt <- model.weights(lm.obj$model)
    wt <- if (is.null(substitute(weights))) rep(1,nrow(data)) else eval(substitute(weights),data,parent.frame())
    ##
    ## XD matrix
    ## lpfunc <- function(delta,fit,dataset,var) {
    ##   dataset[[var]] <- dataset[[var]]+delta
    ##   lpmatrix.lm(fit,dataset)
    ## }
    lpfunc <- if (penalised) function(x,...) {
        newdata <- data
        newdata[[timeVar]] <- x
        predict(lm.obj,newdata,type="lpmatrix")
    } else function(x,fit,data,var) {
        data[[var]] <- x
        lpmatrix.lm(fit,data)
    }
    ## initialise values specific to either delayed entry or interval-censored
    ind0 <- FALSE
    map0 <- 0L
    which0 <- 0
    wt0 <- 0
    ttype <- 0
    transX <- function(X, data) X
    transXD <- function(XD) XD
    smooth <- if (penalised) lm.obj$smooth else NULL
    if (!interval) { # surv.type %in% c("right","counting")
        X <- if (penalised) predict(lm.obj,data,type="lpmatrix") else lpmatrix.lm(lm.obj,data)
        if (link.type=="AH") {
            datat0 <- data
            datat0[[timeVar]] <- 0
            X00 <- if (penalised) predict(lm.obj, datat0, type="lpmatrix")
                   else lpmatrix.lm(lm.obj,datat0)
            index0 <- which.dim(X - X00)
            if (penalised)
                smooth <- lapply(smooth, function(smoothi) {
                    Sindex <- which((1:ncol(X) %in% index0)[smoothi$first.para:smoothi$last.para])
                    para <- range(which((1:ncol(X) %in% smoothi$first.para:smoothi$last.para)[index0]))
                    smoothi$S[[1]] <- smoothi$S[[1]][Sindex,Sindex]
                    smoothi$first.para <- para[1]
                    smoothi$last.para <- para[2]
                    smoothi
                })
            rm(X00)
            transX <- function(X, data) {
                datat0 <- data
                datat0[[timeVar]] <- 0
                Xt0 <- if (penalised) predict(lm.obj,datat0,type="lpmatrix")
                       else lpmatrix.lm(lm.obj,datat0)
                (X - Xt0)[, index0, drop=FALSE]
            }
            transXD <- function(XD) XD[, index0, drop=FALSE]
            init <- init[index0]
        }
        X <- transX(X,data)
        XD <- if (penalised) transXD(grad1(lpfunc,data[[timeVar]],
                                           log.transform=log.time.transform))
              else transXD(grad1(lpfunc,data[[timeVar]],lm.obj,data,timeVar,
                                 log.transform=log.time.transform))
        X1 <- matrix(0,nrow(X),ncol(X))
        X0 <- matrix(0,1,ncol(X))
        if (delayed && all(time0==0)) delayed <- FALSE # CAREFUL HERE: delayed redefined
        if (delayed) {
            ind0 <- time0>0
            map0 <- vector("integer",nrow(X))
            map0[ind0] <- as.integer(1:sum(ind0))
            map0[!ind0] <- NaN
            ##which0 <- which(ind0)
            which0 <- 1:nrow(X)
            which0[!ind0] <- NaN
            data0 <- data[ind0,,drop=FALSE] # data for delayed entry times
            data0[[timeVar]] <- data0[[time0Var]]
            X0 <- if (penalised) transX(predict(lm.obj,data0,type="lpmatrix"), data0)
                  else transX(lpmatrix.lm(lm.obj, data0), data0)
            wt0 <- wt[ind0]
            rm(data0)
        }
    } else { ## interval-censored
        ## ttime <- eventInstance[,1]
        ## ttime2 <- eventInstance[,2]
        ttype <- eventInstance[,3]
        X <- if (penalised) transX(predict(lm.obj,data,type="lpmatrix"), data)
             else transX(lpmatrix.lm(lm.obj,data),data)
        XD <- if (penalised) transXD(grad1(lpfunc,data[[timeVar]],
                                           log.transform=log.time.transform))
              else transXD(grad1(lpfunc,data[[timeVar]],lm.obj,data,timeVar,
                                 log.transform=log.time.transform))
        data0 <- data
        data0[[timeVar]] <- data0[[as.character(time2Expr)]]
        data0[[timeVar]] <- ifelse(data0[[timeVar]]<=0 | data0[[timeVar]]==Inf,NA,data0[[timeVar]])
        X1 <- if (penalised) transX(predict(lm.obj,data0,type="lpmatrix"), data0)
              else transX(lpmatrix.lm(lm.obj, data0), data0)
        X0 <- matrix(0,nrow(X),ncol(X))
        rm(data0)
    }
    if (frailty || copula) {
        Z.formula <- Z
        Z <- model.matrix(Z, data)
        if (ncol(Z)>2) stop("Current implementation only allows for one or two random effects")
        if (ncol(Z)==2) {
            init <- c(init,logtheta1=logtheta,corrtrans=0,logtheta2=logtheta)
        } else init <- c(init, logtheta=unname(logtheta))
    } else {
        Z.formula <- NULL
        Z <- matrix(1,1,1)
    }
    ## check for missing initial values
    if (any(is.na(init) | is.nan(init)))
        stop("Some missing initial values - check that the design matrix is full rank.")
    ## smoothing parameters
    ## cases:
    ##  (1) sp fixed
    ##  (2) sp.init
    ##  (3) use GAM
    no.sp <- is.null(sp)
    if (penalised && no.sp) {
        sp <- if(is.null(lm.obj$full.sp)) lm.obj$sp else lm.obj$full.sp
        if (!is.null(sp.init)) sp <- sp.init
    }
    if (length(control$mle2.control$parscale)==1)
        control$mle2.control$parscale <- rep(control$mle2.control$parscale,length(init))
    if (is.null(names(control$mle2.control$parscale)))
        names(control$mle2.control$parscale) <- names(init)
    if (control$use.gr == FALSE)
        control$optimiser <- "NelderMead"
    args <- list(init=init,X=X,XD=XD,bhazard=bhazard,wt=wt,event=ifelse(event,1,0),time=time,
                 delayed=delayed, interval=interval, X0=X0, wt0=wt0, X1=X1,
                 parscale=control$mle2.control$parscale,
                 kappa=control$kappa.init,outer_optim=control$outer_optim,
                 smooth=if(control$penalty == "logH") smooth else design,
                 sp=sp, reltol_search=control$reltol.search, reltol=control$reltol.final,
                 reltol_outer=control$reltol.outer, trace=control$trace,
                 alpha=alpha, criterion=switch(control$criterion,GCV=1,BIC=2),
                 oldcluster=cluster, frailty=frailty, robust=robust,
                 cluster=if(!is.null(cluster)) as.vector(unclass(factor(cluster))) else NULL,
                 map0 = map0 - 1L, ind0 = ind0, which0 = which0 - 1L, link=link.type, ttype=ttype,
                 RandDist=RandDist, optimiser=control$optimiser, log.time.transform=log.time.transform,
                 type=if (frailty && RandDist=="Gamma") "stpm2_gamma_frailty" else if (frailty && RandDist=="LogN") "stpm2_normal_frailty" else "stpm2",
                 recurrent = recurrent, return_type="optim", transX=transX, transXD=transXD,
                 maxkappa=control$maxkappa, Z=Z, Z.formula = Z.formula, thetaAO = theta.AO,
                 excess=excess, data=data, copula=copula,
                 robust_initial = control$robust_initial, .include=.include,
                 offset=as.vector(offset))
    ## checks on the parameters
    stopifnot(all(dim(args$X) == dim(args$XD)))
    if (!frailty && !copula) {
        stopifnot(length(args$init) == ncol(args$X))
    } else  {
        stopifnot(length(args$init) == ncol(args$X)+ncol(Z))
    }
    if (penalised) args$type <- if (frailty && RandDist=="Gamma") "pstpm2_gamma_frailty"
                                else if (frailty && RandDist=="LogN") "pstpm2_normal_frailty"
                                else "pstpm2"
    if (copula) args$type <- if (penalised) "pstpm2_clayton_copula" else "stpm2_clayton_copula"
    if (frailty || copula) {
        rule <- fastGHQuad::gaussHermiteData(control$nodes)
        args$gauss_x <- rule$x
        args$gauss_w <- rule$w
        args$adaptive <- control$adaptive
        if (ncol(args$Z)>1) {
            control$use.gr <- FALSE
            if(penalised)
                stop("Multiple frailties not implemented for penalised models.")
            args$type <- "stpm2_normal_frailty_2d"
            args$Z <- as.matrix(args$Z)
            ## args$adaptive <- FALSE # adaptive not defined
            ## args$optimiser <- "NelderMead" # gradients not defined
        } else args$Z <- as.vector(args$Z)
    }
    if (penalised) {
        ## penalty function
        pfun <- function(beta,sp) {
            sum(sapply(1:length(lm.obj$smooth),
                       function(i) {
                           smoother <- lm.obj$smooth[[i]]
                           betai <- beta[smoother$first.para:smoother$last.para]
                           sp[i]/2 * betai %*% smoother$S[[1]] %*% betai
                       }))
        }
        negllsp <- function(beta,sp) {
            localargs <- args
            localargs$sp <- sp
            localargs$init <- beta
            localargs$return_type <- "objective"
            negll <- .Call("model_output", localargs, PACKAGE="rstpm2")
            localargs$return_type <- "feasible"
            feasible <- .Call("model_output", localargs, PACKAGE="rstpm2")
            attr(negll,"feasible") <- feasible
            return(negll)
        }
        negll0sp <- function(beta,sp) {
            localargs <- args
            localargs$sp <- sp
            localargs$init <- beta
            localargs$return_type <- "objective0"
            negll <- .Call("model_output", localargs, PACKAGE="rstpm2")
            localargs$return_type <- "feasible"
            feasible <- .Call("model_output", localargs, PACKAGE="rstpm2")
            attr(negll,"feasible") <- feasible
            return(negll)
        }
        ## unused?
        dpfun <- function(beta,sp) {
            deriv <- beta*0
            for (i in 1:length(lm.obj$smooth))
            {
                smoother <- lm.obj$smooth[[i]]
                ind <- smoother$first.para:smoother$last.para
                deriv[ind] <- sp[i] * smoother$S[[1]] %*% beta[ind]
            }
            return(deriv)
        }
        if (control$penalty == "h") {
            ## a current limitation is that the hazard penalty needs to extract the variable names from the smoother objects (e.g. log(time) will not work)
            stopifnot(sapply(lm.obj$smooth,function(obj) obj$term) %in% names(data) ||
                      !is.null(smoother.parameters))
            ## new penalty using the second derivative of the hazard
            design <- smootherDesign(lm.obj,data,smoother.parameters)
            pfun <- function(beta,sp) {
                sum(sapply(1:length(design), function(i) {
                    obj <- design[[i]]
                    s0 <- as.vector(obj$X0 %*% beta)
                    s1 <- as.vector(obj$X1 %*% beta)
                    s2 <- as.vector(obj$X2 %*% beta)
                    s3 <- as.vector(obj$X3 %*% beta)
                    h2 <- (s3+3*s1*s2+s1^3)*exp(s0)
                    sp[i]/2*obj$lambda*sum(obj$w*h2^2)
                }))
            }
            dpfun <- function(beta,sp) {
                if (frailty || copula) beta <- beta[-length(beta)]
                deriv <- beta*0
                for (i in 1:length(design)) {
                    obj <- design[[i]]
                    s0 <- as.vector(obj$X0 %*% beta)
                    s1 <- as.vector(obj$X1 %*% beta)
                    s2 <- as.vector(obj$X2 %*% beta)
                    s3 <- as.vector(obj$X3 %*% beta)
                    h2 <- (s3+3*s1*s2+s1^3)*exp(s0)
                    dh2sq.dbeta <- 2*h2*(exp(s0)*(obj$X3+3*(obj$X1*s2+obj$X2*s1)+3*s1^2*obj$X1)+h2*obj$X0)
                    deriv <- deriv + sp[i]*obj$lambda*colSums(obj$w*dh2sq.dbeta)
                }
                deriv
            }
        }
        gradnegllsp <- function(beta,sp) {
            localargs <- args
            localargs$init <- beta
            localargs$return_type <- "gradient"
            .Call("model_output", localargs, PACKAGE="rstpm2")
        }
        gradnegll0sp <- function(beta,sp) {
            localargs <- args
            localargs$init <- beta
            localargs$return_type <- "gradient0"
            .Call("model_output", localargs, PACKAGE="rstpm2")
        }
        logli <- function(beta) {
            localargs <- args
            localargs$init <- beta
            localargs$return_type <- "li"
            return(.Call("model_output", localargs, PACKAGE="rstpm2"))
        }
        args$logli2 <- function(beta,offset) {
            localargs <- args
            localargs$init <- beta
            localargs$offset <- offset
            localargs$return_type <- "li"
            return(.Call("model_output", localargs, PACKAGE="rstpm2"))
        }
        like <- function(beta) {
            eta <- as.vector(X %*% beta)
            etaD <- as.vector(XD %*% beta)
            h <- link$h(eta,etaD) + bhazard
            H <- link$H(eta)
            ll <- sum(wt[event]*log(h[event])) - sum(wt*H)
            if (delayed) {
                eta0 <- as.vector(X0 %*% beta)
                ## etaD0 <- as.vector(XD0 %*% beta)
                ll <- ll + sum(wt0*link$H(eta0))
            }
            return(ll)
        }
        if (no.sp && !is.null(sp.init)) {
            if(!is.null(lm.obj$full.sp)) lm.obj$sp <- lm.obj$full.sp
            value <- NULL
            while(is.na(value <- negllsp(init,lm.obj$sp)) || !attr(value,"feasible")) {
                lm.call$sp <- lm.obj$sp * 5
                if (no.sp) sp <- lm.call$sp
                ## Unresolved: should we change sp.init if the initial values are not feasible?
                lm.obj <- eval(lm.call)
                if(!is.null(lm.obj$full.sp)) lm.obj$sp <- lm.obj$full.sp
                init <- coef(lm.obj)
                if (frailty || copula)
                    init <- c(init,logtheta=logtheta)
                if (all(lm.obj$sp > 1e5)) break
                ## stop("Initial values not valid and revised sp>1e5")
            }
            args$sp <- lm.obj$sp
        } else args$sp <- sp
        ##     ### Using exterior penalty method for nonlinear constraints: h(t)>=0 or increasing logH(t)
        ##     ### Some initial values should be outside the feasible region
        ##     while(all(XD%*%init>=0)){
        ##       init <- init+0.001
        ##     }
        ##     ### Check initial value
        ##     if(any(XD%*%init<=0)) {
        ##       cat("Some initial values are exactly outside the feasible region of this problem","\n")
        ##     }
        ## MLE
        args$return_type <- if (!no.sp) "optim_fixed"
                            else if (length(sp)>1) "optim_multivariate"
                            else "optim_first"
    } else { # un-penalised models
        negll <- function(beta) {
            localargs <- args
            localargs$return_type <- "objective"
            localargs$init <- beta
            return(.Call("model_output", localargs, PACKAGE="rstpm2"))
        }
        gradnegll <- function(beta) {
            localargs <- args
            localargs$init <- beta
            localargs$return_type <- "gradient"
            return(.Call("model_output", localargs, PACKAGE="rstpm2"))
        }
        fdgradnegll <- function(beta, eps=1e-6) {
            sapply(1:length(beta), function(i) {
                betau <- betal <- beta
                betau[i] <- beta[i]+eps
                betal[i] <- beta[i]-eps
                (negll(betau)-negll(betal))/2/eps
            })
        }
        logli <- function(beta) {
            localargs <- args
            localargs$init <- beta
            localargs$return_type <- "li"
            return(.Call("model_output", localargs, PACKAGE="rstpm2"))
        }
        args$logli2 <- function(beta,offset) {
            localargs <- args
            localargs$init <- beta
            localargs$offset <- offset
            localargs$return_type <- "li"
            return(.Call("model_output", localargs, PACKAGE="rstpm2"))
        }
        parnames(negll) <- parnames(gradnegll) <- names(init)
    }
    ## MLE
    if ((frailty || copula) && !has.init) { # first fit without the frailty
        args2 <- args
        args2$frailty <- args2$copula <- FALSE
        args2$cluster <- NULL
        args2$type <- if (penalised) "pstpm2" else "stpm2"
        localIndex <- 1:(length(args2$init)-1)
        args2$init <- args2$init[localIndex]
        args2$parscale <- args2$parscale[localIndex]
        fit <- .Call("model_output", args2, PACKAGE="rstpm2")
        rm(args2)
        args$init <- c(fit$coef,logtheta)
    }
    fit <- .Call("model_output", args, PACKAGE="rstpm2")
    args$init <- coef <- as.vector(fit$coef)
    if (penalised) {
        args$sp <- sp <- fit$sp <- as.vector(fit$sp)
        edf <- fit$edf
        edf_var<- as.vector(fit$edf_var)
        names(edf_var) <- sapply(lm.obj$smooth,"[[","label")
        negll <- function(beta) negllsp(beta,sp)
        gradnegll <- function(beta) gradnegllsp(beta,sp)
        parnames(negll) <- parnames(gradnegll) <- names(init)
    }
    args$kappa.final <- fit$kappa
    hessian <- fit$hessian
    names(coef) <- rownames(hessian) <- colnames(hessian) <- names(init)
    mle2 <- if (control$use.gr) mle2(negll, coef, vecpar=TRUE, control=control$mle2.control,
                                     gr=gradnegll, ..., eval.only=TRUE)
            else mle2(negll, coef, vecpar=TRUE, control=control$mle2.control, ..., eval.only=TRUE)
    mle2@details$convergence <- if (penalised) 0 else fit$fail # fit$itrmcd
    if (inherits(vcov <- try(solve(hessian,tol=0)), "try-error")) {
        if (control$optimiser=="NelderMead") {
            warning("Non-invertible Hessian")
            mle2@vcov <- matrix(NA,length(coef), length(coef))
        }
        if (control$optimiser!="NelderMead") {
            warning("Non-invertible Hessian - refitting with Nelder-Mead")
            args$optimiser <- "NelderMead"
            fit <- .Call("model_output", args, PACKAGE="rstpm2")
            args$init <- coef <- as.vector(fit$coef)
            args$kappa.final <- fit$kappa
            hessian <- fit$hessian
            names(coef) <- rownames(hessian) <- colnames(hessian) <- names(init)
            mle2 <- mle2(negll, coef, vecpar=TRUE, control=control, ..., eval.only=TRUE)
            mle2@vcov <- if (!inherits(vcov <- try(solve(hessian,tol=0)), "try-error")) vcov else matrix(NA,length(coef), length(coef))
            mle2@details$convergence <- fit$fail # fit$itrmcd
            if (inherits(vcov <- try(solve(hessian,tol=0)), "try-error"))
                warning("Non-invertible Hessian - refitting failed")
        }
    } else {
        mle2@vcov <- vcov
    }
    mle2@details$conv <- mle2@details$convergence
    args$gradnegll = gradnegll
    if (penalised) {
        out <- new("pstpm2",
                   call = mle2@call,
                   call.orig = mle2@call,
                   coef = mle2@coef,
                   fullcoef = mle2@fullcoef,
                   vcov = mle2@vcov,
                   min = mle2@min,
                   details = mle2@details,
                   minuslogl = mle2@minuslogl,
                   method = mle2@method,
                   optimizer = "optim", # mle2@optimizer
                   data = data, # mle2@data, which uses as.list()
                   formula = mle2@formula,
                   xlevels = .getXlevels(mt, mf),
                   ##contrasts = attr(X, "contrasts"),
                   contrasts = NULL, # wrong!
                   logli = logli,
                   ##weights = weights,
                   Call = Call,
                   terms = mt,
                   model.frame = mf,
                   gam = lm.obj,
                   timeVar = timeVar,
                   time0Var = time0Var,
                   timeExpr = timeExpr,
                   time0Expr = time0Expr,
                   like = like,
                   ## fullformula = fullformula,
                   delayed=delayed,
                   frailty = frailty,
                   x = X,
                   xd = XD,
                   termsd = mt, # wrong!
                   y = y,
                   sp = sp,
                   nevent=nevent,
                   link=link,
                   edf=edf,
                   edf_var=edf_var,
                   df=edf,
                   args=args)
        if (robust && !frailty) {
            ## Bread matrix
            bread.mat <- solve(fit$hessian,tol=0)
            ## Meat matirx calculated with individual penalized score functions
            beta.est <- fit$coef
            sp.opt <- fit$sp
            eta <- as.vector(args$X %*% beta.est)
            etaD <- as.vector(XD %*% beta.est)
            h <- link$h(eta,etaD) + bhazard
            H <- link$H(eta)
            gradh <- link$gradh(eta,etaD,args)
            gradH <- link$gradH(eta,args)
            ## right censored data
            score.ind <- t(wt*(gradH - ifelse(event,1/h,0)*gradh)) + dpfun(beta.est, sp.opt)/nrow(gradH)
            meat.mat <- var(t(score.ind))*nrow(gradH)
            out@vcov <- bread.mat %*% meat.mat %*% t(bread.mat)
        }
    } else {
        out <- new("stpm2",
                   call = mle2@call,
                   call.orig = mle2@call,
                   coef = mle2@coef,
                   fullcoef = mle2@fullcoef,
                   vcov = mle2@vcov,
                   min = mle2@min,
                   details = mle2@details,
                   minuslogl = mle2@minuslogl,
                   method = mle2@method,
                   data = as.data.frame(data),
                   formula = mle2@formula,
                   optimizer = "optim",
                   xlevels = .getXlevels(mt, mf),
                   ##contrasts = attr(X, "contrasts"),
                   contrasts = contrasts,
                   logli = logli,
                   ##weights = weights,
                   Call = Call,
                   terms = mt,
                   model.frame = mf,
                   lm = lm.obj,
                   timeVar = timeVar,
                   time0Var = time0Var,
                   timeExpr = timeExpr,
                   time0Expr = time0Expr,
                   delayed = delayed,
                   interval = interval,
                   frailty = frailty,
                   call.formula = formula,
                   x = X,
                   xd = XD,
                   termsd = mt, # wrong!
                   y = y,
                   link=link,
                   args=args)
        if (robust && !frailty) # kludge
          out@vcov <- sandwich.stpm2(out, cluster=cluster)
        attr(out,"nobs") <- nrow(out@x) # for logLik method
    }
    return(out)
}
stpm2 <- function(formula, data, weights=NULL, subset=NULL, coxph.strata=NULL, ...) {
    m <- match.call()
    m[[1L]] <- quote(gsm)
    m$penalised <- FALSE
    out <- eval(m,data,parent.frame())
    out@Call <- match.call()
    out
}
pstpm2 <- function(formula, data, weights=NULL, subset=NULL, coxph.strata=NULL, ...) {
    m <- match.call()
    m[[1L]] <- quote(gsm)
    m$penalised <- TRUE
    out <- eval(m,data,parent.frame())
    out@Call <- match.call()
    out
}

setMethod("update", "stpm2", function(object, ...) {
    object@call = object@Call
    update.default(object, ...)
})
setMethod("show", "stpm2",
          function(object) {
              object@call.orig <- object@Call
              show(as(object,"mle2"))
              })
setClass("summary.stpm2",
         representation(frailty="logical",theta="list",wald="matrix",args="list"),
         contains="summary.mle2")
corrtrans <- function(x) (1-exp(-x)) / (1+exp(-x))
setMethod("summary", "stpm2",
          function(object) {
              newobj <- as(summary(as(object,"mle2")),"summary.stpm2")
              newobj@args <- object@args
              newobj@frailty <- object@frailty
              newobj@call <- object@Call
              if ((object@frailty || object@args$copula) && !is.matrix(object@args$Z)) {
                  coef <- coef(newobj)
                  theta <- exp(coef[nrow(coef),1])
                  se.logtheta <- coef[nrow(coef),2]
                  se.theta <- theta*se.logtheta
                  test.statistic <- (1/se.logtheta)^2
                  p.value <- pchisq(test.statistic,df=1,lower.tail=FALSE)/2
                  newobj@theta <- list(theta=theta, se.theta=se.theta, p.value=p.value)
              } else if (object@frailty && is.matrix(object@args$Z) && ncol(object@args$Z)==2) {
                  coef <- coef(object)
                  index <- (length(coef)-2):length(coef)
                  coef <- coef[index]
                  vcov <- vcov(object)[index,index]
                  rho <- corrtrans(g12 <- coef[2])
                  theta <- c(theta1=exp(coef[1]),corr=rho,theta2=exp(coef[3]))
                  se.theta <- c(theta[1]*sqrt(vcov[1,1]),
                                2*exp(-g12)/(1+exp(-g12))^2*sqrt(vcov[2,2]),
                                theta[3]*sqrt(vcov[3,3]))
                  se.logtheta <- c(theta1=sqrt(vcov[1,1]),
                                   corr=sqrt(vcov[2,2])/coef[2],
                                   theta2=sqrt(vcov[3,3]))
                  test.statistic <- (1/se.logtheta)^2
                  p.value <- pchisq(test.statistic,df=1,lower.tail=FALSE)/2
                  newobj@theta <- list(theta=theta, se.theta=se.theta, p.value=p.value)
              } else newobj@theta <- list()
              newobj@wald <- matrix(NA,1,1) # needed by summary.pstpm2
              newobj })
setMethod("show", "summary.stpm2",
          function(object) {
              show(as(object,"summary.mle2"))
              if (object@frailty || object@args$copula) {
                  if (is.matrix(object@args$Z)) {
                      cat("Random effects model: corr=(1-exp(-corrtrans))/(1+exp(-corrtrans))")
                      cat(sprintf("\ntheta1=%g\tse=%g\tp=%g\n",
                                  object@theta$theta[1],object@theta$se.theta[1],object@theta$p.value[1]))
                      cat(sprintf("\ncorr=%g\tse=%g\tp=%g\n",
                                  object@theta$theta[2],object@theta$se.theta[2],object@theta$p.value[2]))
                      cat(sprintf("\ntheta2=%g\tse=%g\tp=%g\n",
                                  object@theta$theta[3],object@theta$se.theta[3],object@theta$p.value[3]))
                  } else {
                      cat(sprintf("\ntheta=%g\tse=%g\tp=%g\n",
                                  object@theta$theta,object@theta$se.theta,object@theta$p.value))
                  }
              }
          })
setMethod("predictnl", "stpm2",
          function(object,fun,newdata=NULL,link=c("I","log","cloglog","logit"), gd=NULL, ...)
  {
    link <- match.arg(link)
    ## invlinkf <- switch(link,I=I,log=exp,cloglog=cexpexp,logit=expit)
    linkf <- eval(parse(text=link))
    if (is.null(newdata) && !is.null(object@data))
      newdata <- object@data
    localf <- function(coef,...)
      {
        object@fullcoef = coef
        linkf(fun(object,...))
      }
    dm <- numDeltaMethod(object,localf,newdata=newdata,gd=gd,...)
    ## ci <- confint(dm, level = level)
    ## out <- invlinkf(data.frame(Estimate=dm$fit,
    ##                            lower=ci[,1],
    ##                            upper=ci[,2]))
    ## ## cloglog switches the bounds
    ## if (link=="cloglog") 
    ##   out <- data.frame(Estimate=out$Estimate,lower=out$upper,upper=out$lower)
    ## return(out)
    return(dm)
  })
##

predictnl.aft <- function(object,fun,newdata=NULL,link=c("I","log","cloglog","logit"), gd=NULL, ...)
  {
    link <- match.arg(link)
    linkf <- eval(parse(text=link))
    if (is.null(newdata) && !is.null(object@args$data))
      newdata <- object@args$data
    localf <- function(coef,...)
      {
          object@args$init <- object@fullcoef <- coef
        linkf(fun(object,...))
      }
    numDeltaMethod(object,localf,newdata=newdata,gd=gd,...)
  }

setMethod("predictnl", "aft", predictnl.aft)

residuals.stpm2.base <- function(object, type=c("li","gradli")) {
    type <- match.arg(type)
    args <- object@args
    if (type=="li") {
        localargs <- args
        localargs$return_type <- "li"
        out <- as.vector(.Call("model_output", localargs, PACKAGE="rstpm2"))
        names(out) <- rownames(args$X)
    }
    if (type=="gradli") {
        localargs <- args
        localargs$return_type <- "gradli"
        out <- .Call("model_output", localargs, PACKAGE="rstpm2")
        dimnames(out) <- dimnames(args$X)
    }
    return(out)
}    
setMethod("residuals", "stpm2",
          function(object, type=c("li","gradli"))
              residuals.stpm2.base(object=object, type=match.arg(type)))

predict.stpm2.base <- 
          function(object, newdata=NULL,
                   type=c("surv","cumhaz","hazard","density","hr","sdiff","hdiff","loghazard","link","meansurv","meansurvdiff","meanhr","odds","or","margsurv","marghaz","marghr","meanhaz","af","fail","margfail","meanmargsurv","uncured","rmst","probcure","lpmatrix","gradh","gradH","rmstdiff","lpmatrixD"),
                   grid=FALSE, seqLength=300,
                   type.relsurv=c("excess","total","other"), ratetable = survival::survexp.us,
                   rmap, scale=365.24,
                   se.fit=FALSE, link=NULL, exposed=NULL, var=NULL, keep.attributes=FALSE,
                   use.gr=TRUE, level=0.95, n.gauss.quad=100, full=FALSE, ...)
{
    type <- match.arg(type)
    type.relsurv <- match.arg(type.relsurv)
    args <- object@args
    if (type %in% c("fail","margfail")) {
        out <- 1-predict.stpm2.base(object, newdata=newdata,
                                    type=switch(type, fail="surv", margfail="margsurv"),
                                    grid=grid, seqLength=seqLength, type.relsurv=type.relsurv,
                                    ratetable=ratetable, rmap=rmap, scale=scale,
                                    se.fit=se.fit, link=link, exposed=exposed,
                                    var=var, keep.attributes=keep.attributes,
                                    use.gr=use.gr, level=level, n.gauss.quad=n.gauss.quad,
                                    full=full, ...)
        if (se.fit) {temp <- out$lower; out$lower <- out$upper; out$upper <- temp}
        return(out)
    }
    if (is.null(link)) {
        if(args$excess){
            link <- switch(type, surv = "log", cumhaz = "I",
                           hazard = "I", hr = "I", sdiff = "I", hdiff = "I",
                           loghazard = "I", link = "I", odds = "log", or = "log",
                           margsurv = "log", marghaz = "I", marghr = "I",
                           meansurv = "I", meanhr = "I", meanhaz = "I", af = "I",
                           meansurvdiff = "I",
                           fail = "cloglog", uncured = "log", density = "log",
                           rmst = "I", probcure = "cloglog", lpmatrix="I", gradh="I",
                           gradH="I", rmstdiff="I", lpmatrixD="I")
        } else {
            link <- switch(type, surv = "cloglog", cumhaz = "log",
                           hazard = "log", hr = "log", sdiff = "I", hdiff = "I",
                           loghazard = "I", link = "I", odds = "log", or = "log",
                           margsurv = "cloglog", marghaz = "log", marghr = "log",
                           meansurv = "I", meanhr="log", meanhaz = "I", af = "I",
                           meansurvdiff = "I",
                           fail = "cloglog", uncured = "cloglog", density = "log",
                           rmst = "I", probcure = "cloglog", lpmatrix="I", gradh="I",
                           gradH="I", rmstdiff="I", lpmatrixD="I")
        }
    }
    invlinkf <- switch(link,I=I,log=exp,cloglog=cexpexp,logit=expit)
    linkf <- eval(parse(text=link))
    if (type %in% c("uncured","probcure") && is.null(exposed))
        exposed <- function(data) {
            data[[object@timeVar]] <- max(args$time[which(args$event==1)])
            data
        }
    if (is.null(exposed) && is.null(var) & type %in% c("hr","sdiff","hdiff","meansurvdiff","meanhr","or","marghr","af","uncured","probcure","rmstdiff"))
          stop('Either exposed or var required for type in ("hr","sdiff","hdiff","meansurvdiff","meanhr","or","marghr","af","uncured","probcure")')
      if (type %in% c('margsurv','marghaz','marghr','margfail','meanmargsurv') && !object@args$frailty)
          stop("Marginal prediction only for frailty models")
    if (is.null(exposed) && !is.null(var))
        exposed  <- incrVar(var)
    ## exposed is a function that takes newdata and returns the revised newdata
    ## var is a string for a variable that defines a unit change in exposure
    if (is.null(newdata) && type %in% c("hr","sdiff","hdiff","meansurvdiff","meanhr","or","marghr","uncured","rmstdiff"))
        stop("Prediction using type in ('hr','sdiff','hdiff','meansurvdiff','meanhr','or','marghr','uncured','probcure') requires newdata to be specified.")
    calcX <- !is.null(newdata)
    time <- NULL
    if (is.null(newdata)) {
        ##mm <- X <- model.matrix(object) # fails (missing timevar)
        X <- object@x
        XD <- object@xd
        ##y <- model.response(object@model.frame)
        y <- object@y
        time <- as.vector(y[,ncol(y)-1])
        newdata <- as.data.frame(object@data)
    }
    newdata2 <- NULL
    lpfunc <- function(newdata)
        if (inherits(object,"pstpm2")) 
        function(x,...) {
            newdata2 <- newdata
            newdata2[[object@timeVar]] <- x
            predict(object@gam,newdata2,type="lpmatrix")
        } else
            function(x,...) {
                newdata2 <- newdata
                newdata2[[object@timeVar]] <- x
                lpmatrix.lm(object@lm,newdata2)
            }
        ##     function(delta,fit,data,var) {
        ##     data[[var]] <- data[[var]]+delta
        ##     lpmatrix.lm(fit,data)
        ## }
    ## resp <- attr(Terms, "variables")[attr(Terms, "response")] 
    ## similarly for the derivatives
    if (grid) {
      Y <- object@y
      event <- Y[,ncol(Y)]==1 | object@args$interval
      time <- object@data[[object@timeVar]]
      eventTimes <- time[event]
      tt <- seq(min(eventTimes),max(eventTimes),length=seqLength)[-1]
      data.x <- data.frame(tt)
      names(data.x) <- object@timeVar
      newdata[[object@timeVar]] <- NULL
      newdata <- merge(newdata,data.x)
      calcX <- TRUE
    }
    if (type %in% c("rmst","rmstdiff")) {
        stopifnot(object@timeVar %in% names(newdata))
        quad <- gauss.quad(n.gauss.quad)
        a <- 0
        olddata = newdata
        b <- newdata[[object@timeVar]]
        time <- t(outer((b-a)/2,quad$nodes) + (a+b)/2)
        weights <- outer(quad$weights, (b-a)/2)
        rowidx = rep(1:nrow(newdata),each=nrow(time))
        newdata <- newdata[rowidx,,drop=FALSE]
        newdata[[object@timeVar]] <- as.vector(time)
        calcX <- TRUE
    }
    if (args$excess) {
        ## rmap <- substitute(rmap)
      if(type.relsurv != "excess"){
        Sstar <- do.call(survexp, list(substitute(I(timeVar*scale)~1,list(timeVar=as.name(object@timeVar))),
                                       ratetable=ratetable,
                                       scale=scale,
                                       rmap=substitute(rmap),
                                       cohort=FALSE,
                                       data=newdata))
        Hstar <- -log(Sstar)
        if (!is.null(newdata2))
          Sstar2 <- do.call(survexp, list(substitute(I(timeVar*scale)~1,list(timeVar=as.name(object@timeVar))),
                                          ratetable=ratetable,
                                          scale=scale,
                                          rmap=rmap,
                                          cohort=FALSE,
                                          data=newdata2))
        Hstar <- -log(Sstar)
      } 
    }
    if (any(idx <- newdata[object@timeVar] == 0))
        newdata[[object@timeVar]][idx] <- .Machine$double.xmin
    if (calcX)  {
      if (inherits(object, "stpm2")) {
          X <- object@args$transX(lpmatrix.lm(object@lm, newdata), newdata)
          XD <- grad1(lpfunc(newdata),newdata[[object@timeVar]],
                      log.transform=object@args$log.time.transform)
          XD <- object@args$transXD(XD)
      }
      if (inherits(object, "pstpm2")) {
           X <- object@args$transX(predict(object@gam, newdata, type="lpmatrix"), newdata)      
           XD <- object@args$transXD(grad1(lpfunc(newdata),newdata[[object@timeVar]],
                                           log.transform=object@args$log.time.transform))
      }
        ## X <- args$transX(lpmatrix.lm(object@lm, newdata), newdata)
        ## XD <- grad(lpfunc,0,object@lm,newdata,object@timeVar)
        ## XD <- args$transXD(matrix(XD,nrow=nrow(X)))
    }
    ## if (type %in% c("hazard","hr","sdiff","hdiff","loghazard","or","marghaz","marghr","fail","margsurv","meanmargsurv","uncured")) {
    if (is.null(time)) {
        ## how to elegantly extract the time variable?
        ## timeExpr <- 
        ##   lhs(object@call.formula)[[length(lhs(object@call.formula))-1]]
        time <- eval(object@timeExpr,newdata,parent.frame())
        ##
    }
    ## if (any(time == 0)) time[time==0] <- .Machine$double.eps
    ## if (object@delayed && !object@interval) {
    ##   newdata0 <- newdata
    ##   newdata0[[object@timeVar]] <- newdata[[object@time0Var]]
    ##   X0 <- lpmatrix.lm(object@lm, newdata0)
    ##   ## XD0 <- grad(lpfunc,0,object@lm,newdata,object@timeVar)
    ##   ## XD0 <- matrix(XD0,nrow=nrow(X0))
    ## }
    if (type %in% c("hr","sdiff","hdiff","meansurvdiff","meanhr","or","marghr","af","uncured","probcure","rmstdiff")) {
        newdata2 <- exposed(newdata)
      if (inherits(object, "stpm2")) {
          X2 <- object@args$transX(lpmatrix.lm(object@lm, newdata2), newdata2)
          XD2 <- grad1(lpfunc(newdata2),newdata2[[object@timeVar]],
                       log.transform=object@args$log.time.transform)
          XD2 <- object@args$transXD(XD2)
      }
      if (inherits(object, "pstpm2")) {
           X2 <- object@args$transX(predict(object@gam, newdata2, type="lpmatrix"), newdata2)      
           XD2 <- object@args$transXD(grad1(lpfunc(newdata2),newdata2[[object@timeVar]],
                                            log.transform=object@args$log.time.transform))
      }
        ## X2 <- args$transX(lpmatrix.lm(object@lm, newdata2), newdata2)
        ## XD2 <- grad(lpfunc,0,object@lm,newdata2,object@timeVar)
        ## XD2 <- args$transXD(matrix(XD2,nrow=nrow(X)))
    }
    colMeans <- function(x) colSums(x)/apply(x,2,length)
    if (object@frailty && type %in% c("af","meansurvdiff") && args$RandDist=="Gamma" && !object@args$interval && !object@args$delayed) {
        times <- newdata[[object@timeVar]]
        utimes <- sort(unique(times))
        n <- nrow(X)/length(utimes)
        n.cluster <- length(unique(args$cluster))
        newlink <- object@link
        beta <- coef(object)
        npar <- length(beta)
        logtheta <- beta[npar]
        theta <- exp(beta[npar])
        beta <- beta[-npar]
        Hessian <- solve(vcov(object),tol=0)
        eta <- as.vector(X %*% beta)
        eta2 <- as.vector(X2 %*% beta)
        S <- newlink$ilink(eta)
        S2 <- newlink$ilink(eta2)
        H <- -log(S)
        H2 <- -log(S2)
        marg <- function(logtheta,H) (1+exp(logtheta)*H)^(-1/exp(logtheta))
        margS <- marg(logtheta,H)
        margS2 <- marg(logtheta,H2)
        dmarg.dlogtheta <- function(logtheta,H) {
            theta <- exp(logtheta)
            marg(logtheta,H)*(exp(-logtheta)*log(1+theta*H)-H/(1+theta*H))
        }
        ## dmarg.dbeta <- function(logtheta,H,gradH) {
        ##     theta <- exp(logtheta)
        ##     -marg(logtheta,H)*gradH/(1+theta*H)
        ## }
        ## link <- link.PH; eps <- 1e-5; dmarg.dlogtheta(3,.2); (marg(3+eps,.2)-marg(3-eps,.2))/2/eps
        ## eps <- 1e-5; dmarg.dlogtheta(3,.2); (marg(3+eps,.2)-marg(3-eps,.2))/2/eps
        ##
        meanS <- tapply(margS,times,mean)
        meanS2 <- tapply(margS2,times,mean)
        fit <- switch(type,af=1-(1-meanS2)/(1-meanS),meansurvdiff=meanS-meanS2)
        se.fit <- vector("numeric",length(utimes))
        for (i in 1:length(utimes)) {
            index <- which(times==utimes[i])
            newobj <- object
            newobj@args$X <- X[index,,drop=FALSE]
            newobj@args$XD <- XD[index,,drop=FALSE]
            gradli <- residuals(newobj,type="gradli")
            res <- cbind(margS[index]-mean(margS[index]),margS2[index]-mean(margS2[index]))
            res <- apply(res,2,function(col) tapply(col,args$cluster,sum))
            res <- cbind(res, gradli)
            meat <- stats::var(res, na.rm=TRUE)
            colnames(meat) <- rownames(meat) <- c("S","S0", names(beta),"logtheta")
            S.hessian <- cbind(-diag(2)*n/n.cluster,
                               rbind(colSums(margS[index]*(-newlink$gradH(eta[index],list(X=X[index,,drop=FALSE]))/(1+theta*H[index])))/n.cluster,
                                     colSums(margS2[index]*(-newlink$gradH(eta2[index],list(X=X2[index,,drop=FALSE]))/(1+theta*H2[index])))/n.cluster),
                               c(sum(dmarg.dlogtheta(logtheta,H[index]))/n.cluster,
                                 sum(dmarg.dlogtheta(logtheta,H2[index]))/n.cluster))
            par.hessian <- cbind(matrix(0, nrow = npar, ncol = 2), -Hessian / n.cluster)
            bread <- rbind(S.hessian, par.hessian)
            ibread <- solve(bread,tol=0)
            sandwich <- (ibread %*% meat %*% t(ibread) / n.cluster)[1:2, 1:2]
            gradient <- switch(type,
                               af=as.matrix(c( - (1 - meanS2[i]) / (1 - meanS[i]) ^ 2, 1 / (1 - meanS[i])), nrow = 2, ncol = 1),
                               meansurvdiff=matrix(c(1,-1),nrow=2))
            AF.var <- t(gradient) %*% sandwich %*% gradient
            ## S.var <- sandwich[1, 1]
            ## S0.var <- sandwich[2, 2]
            se.fit[i] <- sqrt(AF.var)
        }
        pred <- data.frame(Estimate=fit, lower=fit-1.96*se.fit, upper=fit+1.96*se.fit) 
        if (keep.attributes)
            attr(pred,"newdata") <- newdata
        return(pred)
    }
    if (object@frailty && type %in% c("meanmargsurv") && args$RandDist=="Gamma" && !object@args$interval && !object@args$delayed) {
        times <- newdata[[object@timeVar]]
        utimes <- sort(unique(times))
        n <- nrow(X)/length(utimes)
        n.cluster <- length(unique(args$cluster))
        newlink <- object@link
        beta <- coef(object)
        npar <- length(beta)
        logtheta <- beta[npar]
        theta <- exp(beta[npar])
        beta <- beta[-npar]
        Hessian <- solve(vcov(object),tol=0)
        eta <- as.vector(X %*% beta)
        S <- newlink$ilink(eta)
        H <- -log(S)
        marg <- function(logtheta,H) (1+exp(logtheta)*H)^(-1/exp(logtheta))
        margS <- marg(logtheta,H)
        dmarg.dlogtheta <- function(logtheta,H) {
            theta <- exp(logtheta)
            marg(logtheta,H)*(exp(-logtheta)*log(1+theta*H)-H/(1+theta*H))
        }
        ## eps <- 1e-5; dmarg.dlogtheta(3,.2); (marg(3+eps,.2)-marg(3-eps,.2))/2/eps
        ##
        meanS <- tapply(margS,times,mean)
        fit <- meanS
        se.fit <- vector("numeric",length(utimes))
        for (i in 1:length(utimes)) {
            index <- which(times==utimes[i])
            newobj <- object
            newobj@args$X <- X[index,,drop=FALSE]
            newobj@args$XD <- XD[index,,drop=FALSE]
            ## ## Attempt to use the sandwich estimator for the covariance matrix with the delta method (-> inflated SEs)
            ## res <- residuals(newobj,type="gradli")
            ## meat <- stats::var(res, na.rm=TRUE) # crossprod(res)/n.cluster
            ## colnames(meat) <- rownames(meat) <- c(names(beta),"logtheta")
            ## bread <- -Hessian / n.cluster
            ## ibread <- solve(bread,tol=0)
            ## sandwich <- ibread %*% meat %*% t(ibread) / n.cluster
            ## g <- c(colSums(margS[index]*(-newlink$gradH(eta[index],list(X=X[index,,drop=FALSE]))/(1+theta*H[index])))/n.cluster,
            ##                    sum(dmarg.dlogtheta(logtheta,H[index]))/n.cluster)
            ## se.fit[i] <- sqrt(sum(matrix(g,nrow=1)%*%(sandwich %*% g)))
            gradli <- residuals(newobj,type="gradli")
            res <- tapply(margS[index]-mean(margS[index]),args$cluster,sum)
            res <- cbind(res, gradli)
            meat <- stats::var(res, na.rm=TRUE)
            ## meat <- crossprod(res)/n.cluster
            colnames(meat) <- rownames(meat) <- c("S", names(beta),"logtheta")
            S.hessian <- c(-n/n.cluster,
                               colSums(margS[index]*(-newlink$gradH(eta[index],list(X=X[index,,drop=FALSE]))/(1+theta*H[index])))/n.cluster,
                               sum(dmarg.dlogtheta(logtheta,H[index]))/n.cluster)
            par.hessian <- cbind(matrix(0, nrow = npar, ncol = 1), -Hessian / n.cluster)
            bread <- rbind(S.hessian, par.hessian)
            ibread <- solve(bread,tol=0)
            sandwich <- (ibread %*% meat %*% t(ibread) / n.cluster)[1, 1]
            se.fit[i] <- sqrt(sandwich)
        }
        pred <- data.frame(Estimate=fit, lower=fit-1.96*se.fit, upper=fit+1.96*se.fit) 
        if (keep.attributes)
            attr(pred,"newdata") <- newdata
        return(pred)
    }
    local <-  function (object, newdata=NULL, type="surv", exposed, ...)
    {
        beta <- coef(object)
        tt <- object@terms
        link <- object@link # cf. link for transformation of the predictions
        if (object@frailty || (is.logical(object@args$copula) && object@args$copula)) {
            theta <- exp(beta[length(beta)])
            beta <- beta[-length(beta)]
            if (object@args$RandDist=="LogN") {
                gauss_x <- object@args$gauss_x
                gauss_w <- object@args$gauss_w
                Z <- model.matrix(args$Z.formula, newdata)
                if (ncol(Z)>1) stop("Current implementation only allows for a single random effect")
                Z <- as.vector(Z)
                }
        }
        eta <- as.vector(X %*% beta)
        etaD <- as.vector(XD %*% beta)
        S <- link$ilink(eta)
        h <- link$h(eta,etaD)
        if (!object@args$excess && any(ifelse(is.na(h),FALSE,h<0))) warning(sprintf("Predicted hazards less than zero (n=%i).",sum(ifelse(is.na(h),FALSE,h<0))))
        H = link$H(eta)
        Sigma = vcov(object)
        if (!args$excess) type.relsurv <- "excess" ## ugly hack
        if (type=="link") {
          return(eta)
        }
        if (type=="lpmatrix") {
          return(X)
        }
        if (type=="lpmatrixD") {
          return(XD)
        }
        if (type=="cumhaz") {
            ## if (object@delayed) {
            ##     eta0 <- as.vector(X0 %*% beta)
            ##     etaD0 <- as.vector(XD0 %*% beta)
            ##     H0 <- link$H(eta0)
            ##     return(H - H0)
            ## }
            ## else 
            return(switch(type.relsurv,excess=H,total=H+Hstar,other=Hstar))
        }
        if (type=="density")
            return (S*h)
        if (type=="surv") {
          return(switch(type.relsurv,excess=S,total=S*Sstar,other=Sstar))
        }
        if (type=="fail") {
          return(switch(type.relsurv, excess=1-S,total=1-S*Sstar, other=1-Sstar))
        }
        if (type=="odds") { # delayed entry?
          return((1-S)/S)
        }
        if (type=="sdiff")
          return(link$ilink(as.vector(X2 %*% beta)) - S)
        if (type=="hazard") {
          return(h)
        }
        if (type=="loghazard") {
            return(log(h))
        }
        if (type=="hdiff") {
            eta2 <- as.vector(X2 %*% beta)
            etaD2 <- as.vector(XD2 %*% beta)
            h2 <- link$h(eta2,etaD2)
            return(h2 - h)
        }
        if (type=="uncured") {
            S2 <- link$ilink(as.vector(X2 %*% beta))
            return((S-S2)/(1-S2))
        }
        if (type=="probcure") {
            S2 <- link$ilink(as.vector(X2 %*% beta))
            return(S2/S)
        }
        if (type=="hr") {
            eta2 <- as.vector(X2 %*% beta)
            etaD2 <- as.vector(XD2 %*% beta)
            h2 <- link$h(eta2,etaD2)
            return(h2/h)
        }
        if (type=="or") {
            S2 <- link$ilink(as.vector(X2 %*% beta)) 
            return((1-S2)/S2/((1-S)/S))
        }
        if (type=="meansurv") {
            return(tapply(S,newdata[[object@timeVar]],mean))
        }
        if (type=="meanhaz") {
            return(tapply(S*h,newdata[[object@timeVar]],sum)/tapply(S,newdata[[object@timeVar]],sum))
        }
        if (type=="meanhr") {
            eta2 <- as.vector(X2 %*% beta)
            etaD2 <- as.vector(XD2 %*% beta)
            h2 <- link$h(eta2,etaD2)
            S2 <- link$ilink(eta2)
            return(tapply(S2*h2,newdata[[object@timeVar]],sum)/tapply(S2,newdata[[object@timeVar]],sum) / (tapply(S*h,newdata[[object@timeVar]],sum)/tapply(S,newdata[[object@timeVar]],sum)))
        }
        if (type=="meansurvdiff") {
            eta2 <- as.vector(X2 %*% beta)
            S2 <- link$ilink(eta2)
            return(tapply(S2,newdata[[object@timeVar]],mean) - tapply(S,newdata[[object@timeVar]],mean))
        }
        if (type=="af") {
            eta2 <- as.vector(X2 %*% beta)
            S2 <- link$ilink(eta2)
            meanS <- tapply(S,newdata[[object@timeVar]],mean)
            meanS2 <- tapply(S2,newdata[[object@timeVar]],mean)
            if (object@frailty) {
                if (object@args$RandDist=="Gamma") {
                meanS <- tapply((1+theta*(-log(S)))^(-1/theta), newdata[[object@timeVar]], mean)
                meanS2 <- tapply((1+theta*(-log(S2)))^(-1/theta), newdata[[object@timeVar]], mean)
                } else {
                    meanS <- tapply(sapply(1:length(gauss_x),
                                           function(i) link$ilink(eta+Z*sqrt(2)*sqrt(theta)*gauss_x[i])) %*%
                                    gauss_w / sqrt(pi),
                                    newdata[[object@timeVar]],
                                    mean)
                    meanS2 <- tapply(sapply(1:length(gauss_x),
                                           function(i) link$ilink(eta2+Z*sqrt(2)*sqrt(theta)*gauss_x[i])) %*%
                                    gauss_w / sqrt(pi),
                                    newdata[[object@timeVar]],
                                    mean)
                    }
                }
            return((meanS2 - meanS)/(1-meanS))
        }
        if (type=="meanmargsurv") {
            stopifnot(object@frailty && object@args$RandDist %in% c("Gamma","LogN"))
            if (object@args$RandDist=="Gamma")
                return(tapply((1+theta*H)^(-1/theta), newdata[[object@timeVar]], mean))
            if (object@args$RandDist=="LogN") {
                return(tapply(sapply(1:length(gauss_x),
                                     function(i) link$ilink(eta+Z*sqrt(2)*sqrt(theta)*gauss_x[i])) %*%
                              gauss_w / sqrt(pi),
                              newdata[[object@timeVar]],
                              mean))
            }
        }
        if (type=="margsurv") {
            stopifnot(object@args$frailty && object@args$RandDist %in% c("Gamma","LogN"))
            if (object@args$RandDist=="Gamma")
                return((1+theta*H)^(-1/theta))
            if (object@args$RandDist=="LogN") {
                return(sapply(1:length(gauss_x),
                              function(i) link$ilink(eta+Z*sqrt(2)*sqrt(theta)*gauss_x[i])) %*%
                       gauss_w / sqrt(pi))
            }
        }
        if (type=="marghaz") {
            stopifnot(object@frailty && object@args$RandDist %in% c("Gamma","LogN"))
            if (object@args$RandDist=="Gamma") {
                ## margsurv <- (1+theta*H)^(-1/theta)
                ## return(h*margsurv^theta)
                return(h/(1+H*theta))
            }
            if (object@args$RandDist=="LogN") {
                return(sapply(1:length(gauss_x),
                              function(i) link$h(eta+Z*sqrt(2)*sqrt(theta)*gauss_x[i],etaD)) %*%
                       gauss_w / sqrt(pi))
            }
        }
        if (type=="marghr") {
            stopifnot(object@frailty && object@args$RandDist %in% c("Gamma","LogN"))
            eta2 <- as.vector(X2 %*% beta)
            etaD2 <- as.vector(XD2 %*% beta)
            if (object@args$RandDist=="Gamma") {
                H2 <- link$H(eta2)
                h2 <- link$h(eta2,etaD2)
                margsurv <- (1+theta*H)^(-1/theta)
                marghaz <- h*margsurv^theta
                margsurv2 <- (1+theta*H2)^(-1/theta)
                marghaz2 <- h2*margsurv2^theta
            }
            if (object@args$RandDist=="LogN") {
                marghaz <- sapply(1:length(gauss_x),
                                  function(i) as.vector(link$h(eta+Z*sqrt(2)*sqrt(theta)*gauss_x[i],etaD))) %*%
                                  gauss_w / sqrt(pi)
                marghaz2 <- sapply(1:length(gauss_x),
                                   function(i) as.vector(link$h(eta2+Z*sqrt(2)*sqrt(theta)*gauss_x[i],etaD2))) %*%
                                       gauss_w / sqrt(pi)
            }
            return(marghaz2/marghaz)
        }
        if (type=="rmst") {
            return(tapply(S*weights,rowidx,sum))
        }
        if (type=="rmstdiff") {
            S2 <- link$ilink(as.vector(X2 %*% beta))
            return(tapply(S*weights,rowidx,sum) - tapply(S2*weights,rowidx,sum))
        }
        if (type=="gradh")
            return(link$gradh(eta,etaD,list(X=X,XD=XD)))
        if (type=="gradH")
            return(link$gradH(eta,list(X=X)))
    }
    if (!se.fit) {
        out <- local(object,newdata,type=type,exposed=exposed,  ...)
    } else {
      gd <- NULL
      beta <- coef(object)
      ## calculate gradients for some of the estimators
      if (use.gr) {
          colMeans <- function(x) apply(x,2,mean)
          collapse <- function(gd)
              do.call("cbind",tapply(1:nrow(gd), newdata[[object@timeVar]], function(index) colMeans(gd[index, ,drop=FALSE])))
          collapse1 <- function(S)
              as.vector(tapply(S, newdata[[object@timeVar]], mean))
          fd <- function(f,x,eps=1e-5)
              t(sapply(1:length(x),
                       function(i) {
                           upper <- lower <- x
                           upper[i]=x[i]+eps
                           lower[i]=x[i]-eps
                           (f(upper)-f(lower))/2/eps
                       }))
          div <- function(mat,v) t(t(mat)/v)
          if (type=="hazard" && link %in% c("I","log")) {
              ## Case: frailty model (assumes baseline hazard for frailty=1)
              betastar <- if(args$frailty || args$copula) beta[-length(beta)] else beta
              gd <- switch(link,
                           I=t(object@link$gradh(X %*% betastar, XD %*% betastar, list(X=X, XD=XD))),
                           log=t(object@link$gradh(X %*% betastar, XD %*% betastar, list(X=X, XD=XD))/
                               object@link$h(X %*% betastar, XD %*% betastar)))
          }
          if (type=="meanhaz" && !object@frailty && link %in% c("I","log")) {
              betastar <- if(args$frailty || args$copula) beta[-length(beta)] else beta
              h <- as.vector(object@link$h(X %*% betastar, XD %*% betastar))
              S <- as.vector(object@link$ilink(X %*% betastar))
              gradS <- object@link$gradS(X%*% betastar,X)
              gradh <- object@link$gradh(X %*% betastar,XD %*% betastar,list(X=X,XD=XD))
              gd <- if(link=="I") collapse(gradh*S + gradS*h) else
                                                                  div(collapse(gradh*S + h*gradS),collapse1(h*S))-div(collapse(gradS),collapse1(S))
          }
          if (type=="meanhr" && !object@frailty && link == "log") {
              betastar <- if(args$frailty|| args$copula) beta[-length(beta)] else beta
              h <- as.vector(object@link$h(X %*% betastar, XD %*% betastar))
              S <- as.vector(object@link$ilink(X %*% betastar))
              gradS <- object@link$gradS(X %*% betastar,X)
              gradh <- object@link$gradh(X %*% betastar,XD %*% betastar,list(X=X,XD=XD))
              h2 <- as.vector(object@link$h(X2 %*% betastar, XD2 %*% betastar))
              S2 <- as.vector(object@link$ilink(X2 %*% betastar))
              gradS2 <- object@link$gradS(X2 %*% betastar,X2)
              gradh2 <- object@link$gradh(X2 %*% betastar,XD2 %*% betastar,list(X=X2,XD=XD2))
              gd <- div(collapse(gradh2*S2 + h2*gradS2),collapse1(h2*S2)) - div(collapse(gradS2),collapse1(S2)) - div(collapse(gradh*S + h*gradS),collapse1(h*S)) + div(collapse(gradS),collapse1(S))
              ## fd(function(beta) {fit=object; fit@fullcoef = beta; log(predict(fit,type="meanhr",newdata=newdata, var=var, exposed=exposed))},betastar)
          }
          if (type=="meansurv" && !object@frailty && link=="I") {
              gd <- collapse(object@link$gradS(X%*% beta,X))
          }
          if (type=="meansurvdiff" && !object@frailty) {
              gd <- collapse(object@link$gradS(X2%*% beta,X2) - object@link$gradS(X%*% beta,X))
          }
          if (type=="margsurv" && link %in% c("I","cloglog") && args$RandDist=="Gamma") {
              theta <- exp(beta[length(beta)])
              betastar <- beta[-length(beta)]
              eta <- as.vector(X %*% betastar)
              H <- as.vector(object@link$H(eta))
              gradH <- object@link$gradH(eta,list(X=X))
              S0 <- 1+theta*H
              margS <- S0^(-1/theta)
              ## This depends on the transformation link
              if (link=="I")
                  gd <- t(cbind(-margS * gradH/(1+theta*H),
                                margS*(1/theta*log(1+theta*H)-H/(1+theta*H))))
              if (link=="cloglog")
                  gd <- t(cbind(-(theta^2*S0^(-theta-1)*gradH/(S0^(-theta)*(-theta)*log(S0))),
                                (theta*log(S0)*S0^(-theta)-theta^2*H*S0^(-1-theta))/(S0^(-theta)*(-theta*log(S0)))))
          }
          if (type=="marghaz" && link %in% c("I","log") && args$RandDist=="Gamma") {
              theta <- exp(beta[length(beta)])
              betastar <- beta[-length(beta)]
              eta <- as.vector(X %*% betastar)
              etaD <- as.vector(XD %*% betastar)
              H <- as.vector(object@link$H(eta))
              h <- as.vector(object@link$h(eta,etaD))
              gradH <- object@link$gradH(eta,list(X=X))
              gradh <- object@link$gradh(eta,etaD,list(X=X,XD=XD))
              S0 <- 1+theta*H
              margS <- S0^(-1/theta)
              ## This depends on the transformation link
              if (link=="I")
                  gd <- t(cbind((S0*gradh-theta*h*gradH)/(theta^2*H^2+S0),
                                -(theta*H*h/theta^2*H^2+S0)))
              if (link=="log")
                  gd <- t(cbind((S0*gradh-theta*h*gradH)/(S0*h),
                                -theta*H/S0))
          }
          if (type=="af" && !object@frailty) {
              meanS <- collapse1(as.vector(object@link$ilink(X%*%beta)))
              meanS2 <- collapse1(as.vector(object@link$ilink(X2%*%beta)))
              gradS <- collapse(object@link$gradS(X%*%beta,X))
              gradS2 <- collapse(object@link$gradS(X2%*%beta,X2))
              gd <- t((t(gradS2-gradS)*(1-meanS) + (meanS2-meanS)*t(gradS))/(1-meanS)^2)
              ## check
              ## fd(function(beta) collapse1(as.vector(object@link$ilink(X %*% beta))), beta) - gradS # ok
              ## fd(function(beta) collapse1(as.vector(object@link$ilink(X2 %*% beta))), beta) - gradS2 # ok
              ## fd(function(beta) collapse1(as.vector(object@link$ilink(X2 %*% beta)-object@link$ilink(X %*% beta)))/collapse1(1-as.vector(object@link$ilink(X %*% beta))),beta) - gd # ok
          }
      }
      pred <- predictnl(object,local,link=link,newdata=newdata,type=type,gd=if (use.gr) gd else NULL,
                exposed=exposed,...) 
      ci <- confint.predictnl(pred, level = level)
      out <- data.frame(Estimate=pred$fit,
                        lower=ci[,1],
                        upper=ci[,2])
      if (link=="cloglog") 
          out <- data.frame(Estimate=out$Estimate,lower=out$upper,upper=out$lower)
      out <- invlinkf(out)
    }
    if (keep.attributes || full) {
      if (type %in% c("rmst","rmstdiff"))
          newdata = olddata
      if (type %in% c("hr","sdiff","hdiff","meansurvdiff","meanhr","or","marghr","af","rmstdiff"))
          newdata <- exposed(newdata)
      if (type %in% c("meansurv","meanhr","meansurvdiff","meanhaz","meanmargsurv","af")) {
          newdata <- data.frame(time=unique(newdata[[object@timeVar]]))
          names(newdata) <- object@timeVar
      }
      if (full) {
          if (inherits(out, "AsIs"))
              class(out) <- setdiff(class(out),"AsIs")
          out <- if(is.data.frame(out)) cbind(newdata,out) else cbind(newdata, data.frame(Estimate=out))
      }
      else attr(out,"newdata") <- newdata
    }
    return(out)
}

predict.cumhaz <-
          function(object, newdata=NULL)
{
    args <- object@args
    lpfunc <- function(newdata)
        if (inherits(object,"pstpm2"))
        function(x,...) {
            newdata2 <- newdata
            newdata2[[object@timeVar]] <- x
            predict(object@gam,newdata2,type="lpmatrix")
        } else
            function(x,...) {
                newdata2 <- newdata
                newdata2[[object@timeVar]] <- x
                lpmatrix.lm(object@lm,newdata2)
            }
    if (inherits(object, "stpm2")) {
          X <- object@args$transX(lpmatrix.lm(object@lm, newdata), newdata)
      }
    if (inherits(object, "pstpm2")) {
           X <- object@args$transX(predict(object@gam, newdata, type="lpmatrix"), newdata)
      }
    link <- object@link # cf. link for transformation of the predictions
    eta <- as.vector(X %*% beta)
    link$H(eta)
}

setMethod("predict", "stpm2",
          function(object,newdata=NULL,
                   type=c("surv","cumhaz","hazard","density","hr","sdiff","hdiff","loghazard","link","meansurv","meansurvdiff","meanhr","odds","or","margsurv","marghaz","marghr","meanhaz","af","fail","margfail","meanmargsurv","uncured","rmst","probcure","lpmatrix","gradh","gradH","rmstdiff","lpmatrixD"),
                   grid=FALSE,seqLength=300,
                   type.relsurv=c("excess","total","other"), scale=365.24, rmap, ratetable=survival::survexp.us,
                   se.fit=FALSE,link=NULL,exposed=NULL,var=NULL,keep.attributes=FALSE,use.gr=TRUE,level=0.95,n.gauss.quad=100,full=FALSE,...) {
              type <- match.arg(type)
              type.relsurv <- match.arg(type.relsurv)
              predict.stpm2.base(object, newdata=newdata, type=type, grid=grid, seqLength=seqLength, se.fit=se.fit,
                                 link=link, exposed=exposed, var=var, keep.attributes=keep.attributes, use.gr=use.gr,level=level,
                                 type.relsurv=type.relsurv, scale=scale, rmap=rmap, ratetable=ratetable, n.gauss.quad=n.gauss.quad, full=full, ...)
          })

##`%c%` <- function(f,g) function(...) g(f(...)) # function composition
plot.meansurv <- function(x, y=NULL, times=NULL, newdata=NULL, type="meansurv", exposed=NULL, var=NULL, add=FALSE, ci=!add, rug=!add, recent=FALSE,
                          xlab=NULL, ylab=NULL, lty=1, line.col=1, ci.col="grey", seqLength=301, ...) {
    ## if (is.null(times)) stop("plot.meansurv: times argument should be specified")
    if (is.null(newdata)) newdata <- as.data.frame(x@data)
    if (is.null(times)) {
      Y <- x@y
      event <- Y[,ncol(Y)]==1 | x@args$interval
      time <- x@data[[x@timeVar]]
      eventTimes <- time[event]
      times <- seq(min(eventTimes),max(eventTimes),length=seqLength)[-1]
      ## data.x <- data.frame(X)
      ## names(data.x) <- object@timeVar
      ## newdata <- merge(newdata,data.x)
    }
    times <- times[times !=0]
    if (recent) {
        newdata <- do.call("rbind",
                           lapply(times, 
                                  function(time) {
                                      newd <- newdata
                                      newd[[x@timeVar]] <- newdata[[x@timeVar]]*0+time
                                      newd
                                  }))
        pred <- predict(x, newdata=newdata, type=type, keep.attributes=TRUE, se.fit=ci, exposed=exposed, var=var) # requires recent version
        if (type=="meansurv") {
          pred <- if (ci) rbind(c(Estimate=1,lower=1,upper=1),pred) else c(1,pred)
          times <- c(0,times)
        }
        if (type=="meansurvdiff") {
          pred <- if (ci) rbind(c(Estimate=0,lower=0,upper=0),pred) else c(0,pred)
          times <- c(0,times)
        }
    } else {
        pred <- lapply(times, 
                       function(time) {
                           newdata[[x@timeVar]] <- newdata[[x@timeVar]]*0+time
                           predict(x, newdata=newdata, type=type, se.fit=ci, grid=FALSE, exposed=exposed, var=var, keep.attributes=TRUE)
                       })
        pred <- do.call("rbind", pred)
        if (type=="meansurv")  {
            pred <- if (ci) rbind(c(Estimate=1,lower=1,upper=1),pred) else c(1,unlist(pred))
            times <- c(0,times)
            }
        if (type=="meansurvdiff")  {
            pred <- if (ci) rbind(c(Estimate=0,lower=0,upper=0),pred) else c(0,unlist(pred))
            times <- c(0,times)
            }
        }
    if (is.null(xlab)) xlab <- deparse(x@timeExpr)
    if (is.null(ylab))
        ylab <- switch(type,
                       meansurv="Mean survival",
                       af="Attributable fraction",
                       meansurvdiff="Difference in mean survival",
                       meanhaz="Mean hazard",
                       meanhr="Mean hazard ratio")
    if (!add) matplot(times, pred, type="n", xlab=xlab, ylab=ylab, ...)
    if (ci) {
        polygon(c(times,rev(times)),c(pred$lower,rev(pred$upper)),col=ci.col,border=ci.col)
        lines(times,pred$Estimate,col=line.col,lty=lty,...)
    } else {
        lines(times,pred,col=line.col,lty=lty,...)
    }
    if (rug) {
        Y <- x@y
        eventTimes <- Y[Y[,ncol(Y)]==1,ncol(Y)-1]
        rug(eventTimes,col=line.col)
    }
    return(invisible(y))
}

plot.stpm2.base <- 
          function(x,y,newdata=NULL,type="surv",
                   xlab=NULL,ylab=NULL,line.col=1,ci.col="grey",lty=par("lty"),log="",
                   add=FALSE,ci=!add,rug=!add,
                   var=NULL,exposed=NULL,times=NULL,
                   type.relsurv=c("excess","total","other"), ratetable = survival::survexp.us, rmap, scale=365.24, ...) {
              if (type %in% c("meansurv","meansurvdiff","af","meanhaz","meanhr")) {
                  return(plot.meansurv(x,times=times,newdata=newdata,type=type,xlab=xlab,ylab=ylab,line.col=line.col,ci.col=ci.col,
                                       lty=lty,add=add,ci=ci,rug=rug, exposed=exposed, var=var, ...))
              }
              if (is.null(newdata)) stop("newdata argument needs to be specified")
              y <- predict(x,newdata,type=switch(type,fail="surv",margfail="margsurv",type),var=var,exposed=exposed,
                           grid=!(x@timeVar %in% names(newdata)), se.fit=ci, keep.attributes=TRUE, 
                           type.relsurv=type.relsurv, ratetable=ratetable, rmap=rmap, scale=scale)
              if (type %in% c("fail","margfail")) {
                  if (ci) {
                      y$Estimate <- 1-y$Estimate
                      lower <- y$lower
                      y$lower=1-y$upper
                      y$upper=1-lower
                      } else y <- structure(1-y,newdata=attr(y,"newdata"))
              }
              if (is.null(xlab)) xlab <- deparse(x@timeExpr)
              if (is.null(ylab))
                  ylab <- switch(type,hr="Hazard ratio",hazard="Hazard",surv="Survival",density="Density",
                                 sdiff="Survival difference",hdiff="Hazard difference",cumhaz="Cumulative hazard",
                                 loghazard="log(hazard)",link="Linear predictor",meansurv="Mean survival",
                                 meansurvdiff="Difference in mean survival",meanhr="Mean hazard ratio",
                                 odds="Odds",or="Odds ratio",
                                 margsurv="Marginal survival",marghaz="Marginal hazard",marghr="Marginal hazard ratio", haz="Hazard",fail="Failure",
                                 meanhaz="Mean hazard",margfail="Marginal failure",af="Attributable fraction",meanmargsurv="Mean marginal survival",
                                 uncured="Uncured distribution",
                                 rmst="Restricted mean survival time",
                                 rmstdiff="Restricted mean survival time difference")
              xx <- attr(y,"newdata")
              xx <- eval(x@timeExpr,xx) # xx[,ncol(xx)]
              ## remove NaN
              if (any(is.nan(as.matrix(y)))) {
                  idx <- is.nan(y[,1])
                  for (j in 2:ncol(y))
                      idx <- idx | is.nan(y[,j])
                  y <- y[!idx,,drop=FALSE]
                  xx <- xx[!idx]
              }
              if (!add) matplot(xx, y, type="n", xlab=xlab, ylab=ylab, log=log, ...)
              if (ci) {
                  polygon(c(xx,rev(xx)), c(y[,2],rev(y[,3])), col=ci.col, border=ci.col)
                  lines(xx,y[,1],col=line.col,lty=lty,...)
              } else lines(xx,y,col=line.col,lty=lty,...)
              if (rug) {
                  Y <- x@y
                  eventTimes <- Y[Y[,ncol(Y)]==1,ncol(Y)-1]
                  rug(eventTimes,col=line.col)
              }
              return(invisible(y))
          }
setMethod("plot", signature(x="stpm2", y="missing"),
          function(x,y,newdata=NULL,type="surv",
                   xlab=NULL,ylab=NULL,line.col=1,ci.col="grey",lty=par("lty"),
                   add=FALSE,ci=!add,rug=!add,
                   var=NULL,exposed=NULL,times=NULL,
                   type.relsurv=c("excess","total","other"), ratetable = survival::survexp.us, 
                   rmap, scale=365.24,...)
              plot.stpm2.base(x=x, y=y, newdata=newdata, type=type, xlab=xlab,
                              ylab=ylab, line.col=line.col, ci.col=ci.col, lty=lty, add=add,
                              ci=ci, rug=rug, var=var, exposed=exposed, times=times, 
                              type.relsurv=type.relsurv, ratetable=ratetable, rmap=rmap, scale=scale,...)
          )
lines.stpm2 <- 
          function(x,newdata=NULL,type="surv",
                   col=1,ci.col="grey",lty=par("lty"),
                   ci=FALSE,rug=FALSE,
                   var=NULL,exposed=NULL,times=NULL,
                   type.relsurv=c("excess","total","other"), ratetable = survival::survexp.us, 
                   rmap, scale=365.24, ...)
              plot.stpm2.base(x=x, newdata=newdata, type=type, 
                              line.col=col, ci.col=ci.col, lty=lty, add=TRUE,
                              ci=ci, rug=rug, var=var, exposed=exposed, times=times, 
                              type.relsurv=type.relsurv, ratetable=ratetable, rmap=rmap, scale=scale,...)
setMethod("lines", signature(x="stpm2"), lines.stpm2)
eform <- function (object, ...) 
  UseMethod("eform")
setGeneric("eform")
eform.stpm2 <- function (object, parm, level = 0.95, method = c("Profile","Delta"), 
                    name = "exp(beta)", ...) 
          {
              method <- match.arg(method)
              if (missing(parm)) 
                  parm <- TRUE
              if (object@args$robust) {
                  ## Profile likelihood not defined for robust variance: using delta method
                  method <- "Delta"
              }
              estfun <- switch(method, Profile = confint, Delta = stats::confint.default)
              val <- exp(cbind(coef = coef(object), estfun(object, level = level)))
              colnames(val) <- c(name, colnames(val)[-1])
              val[parm, ]
          }
setMethod("eform", signature(object="stpm2"), eform.stpm2)
eform.default <- function(object, parm, level = 0.95, method=c("Delta","Profile"), name="exp(beta)", ...) {
  method <- match.arg(method)
  if (missing(parm))
      parm <- TRUE
  if (method == "Profile") class(object) <- c(class(object),"glm")
  estfun <- switch(method, Profile = confint, Delta = stats::confint.default)
  val <- exp(cbind(coef = coef(object), estfun(object, level = level)))
  colnames(val) <- c(name,colnames(val)[-1])
  val[parm, ]
}

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
smootherDesign <- function(gamobj,data,parameters = NULL) {
    d <- data[1,,drop=FALSE] ## how to get mean prediction values, particularly for factors?
    makepred <- function(var,inverse) {
        function(value) {
            d <- d[rep(1,length(value)),]
            d[[var]] <- inverse(value)
            predict(gamobj,newdata=d,type="lpmatrix")
        }
    }
    smoother.names <- sapply(gamobj$smooth, function(obj) obj$term)
    lapply(1:length(gamobj$smooth), function(i) {
        smoother <- gamobj$smooth[[i]]
        if (is.null(parameters)) {
            var <- smoother$term
            stopifnot(var %in% names(data))
            transform <- I
            inverse <- I
        } else {
            j <- match(smoother$term,names(parameters))
            stopifnot(!is.na(j))
            var <- parameters[[j]]$var
            transform <- parameters[[j]]$transform
            inverse <- parameters[[j]]$inverse
        }
        pred <- makepred(var,inverse)
        derivativeDesign(pred,
                         lower=transform(min(data[[var]])),
                         upper=transform(max(data[[var]])))
    })
}
## TODO: If we transform a smoother (e.g. log(time)), we can use information on
## (i) the variable name, (ii) the transform and (iii) the inverse transform.



## penalised stpm2
setOldClass("gam")
setClass("pstpm2", representation(xlevels="list",
                                  contrasts="listOrNULL",
                                  terms="terms",
                                  logli="function",
                                  gam="gam",
                                  timeVar="character",
                                  time0Var="character",
                                  timeExpr="nameOrcall",
                                  time0Expr="nameOrcallOrNULL",
                                  like="function",
	                          model.frame="list",
	                          ## fullformula="formula",
                                  delayed="logical",
                                  frailty="logical",
                                  x="matrix",
                                  xd="matrix",
                                  termsd="terms",
                                  Call="call",
                                  y="Surv",
                                  sp="numeric",
                                  nevent="numeric",
                                  link="list",
                                  edf="numeric",
                                  edf_var="numeric",
                                  df="numeric",
                                  args="list"),
         contains="mle2")

## Could this inherit from summary.stpm2?
setClass("summary.pstpm2", representation(pstpm2="pstpm2",frailty="logical",theta="list",wald="matrix"), contains="summary.mle2")
setMethod("summary", "pstpm2",
          function(object) {
              newobj <- as(summary(as(object,"mle2")),"summary.pstpm2")
              newobj@pstpm2 <- object
              newobj@frailty <- object@frailty
              newobj@call <- object@Call
              if (object@frailty) {
                  coef <- coef(newobj)
                  theta <- exp(coef[nrow(coef),1])
                  se.logtheta <- coef[nrow(coef),2]
                  se.theta <- theta*se.logtheta
                  test.statistic <- (1/se.logtheta)^2
                  p.value <- pchisq(test.statistic,df=1,lower.tail=FALSE)/2
                  newobj@theta <- list(theta=theta, se.theta=se.theta, p.value=p.value)
              } else newobj@theta <- list()
              vcov1 <- vcov(object)
              coef1 <- coef(object)
              ## Wald test for the smoothers
              wald <- t(sapply(names(object@edf_var), function(name) {
                  i <- grep(name,colnames(vcov1),fixed=TRUE)
                  statistic <- as.vector(coef1[i] %*% solve(vcov1[i,i],tol=0) %*% coef1[i])
                  edf <- object@edf_var[name]
                  c(statistic=statistic,ncoef=length(i),edf=edf,p.value=pchisq(statistic, edf, lower.tail=FALSE))
              }))
              colnames(wald) <- c("Wald statistic","Number of coef","Effective df","P value")
              newobj@wald <- wald
              newobj })
setMethod("show", "summary.pstpm2",
          function(object) {
              show(as(object,"summary.mle2"))
              cat(sprintf("\nEffective df=%g\n",object@pstpm2@edf))
              printCoefmat(object@wald)
              if (object@frailty)
                  cat(sprintf("\ntheta=%g\tse=%g\tp=%g\n",
                              object@theta$theta,object@theta$se.theta,object@theta$p.value))
          })
setMethod("AICc", "pstpm2",
          function (object, ..., nobs=NULL, k=2)  {
              L <- list(...)
              if (length(L)) {
                  L <- c(list(object),L)
                  if (is.null(nobs)) {
                      nobs <- sapply(L,nobs)
                  }
                  if (length(unique(nobs))>1)
                      stop("nobs different: must have identical data for all objects")
                  val <- sapply(L, AICc, nobs=nobs, k=k)
                  df <- sapply(L,attr,"edf")
                  data.frame(AICc=val,df=df)
              } else {
                  df <- attr(object,"edf")
                  if (is.null(nobs)) nobs <- object@nevent
                  c(-2*logLik(object)+k*df+k*df*(df+1)/(nobs-df-1))
              }
          })
setMethod("qAICc", "pstpm2",
          function (object, ..., nobs = NULL, dispersion = 1, k = 2)  {
              L <- list(...)
              if (length(L)) {
                  L <- c(list(object),L)
                  if (is.null(nobs)) {
                      nobs <- sapply(L,nobs)
                  }
                  if (length(unique(nobs))>1)
                      stop("nobs different: must have identical data for all objects")
                  val <- sapply(L, qAICc, nobs=nobs,dispersion=dispersion,k=k)
                  df <- sapply(L,attr,"edf")
                  data.frame(qAICc=val,df=df)
              } else {
                  df <- attr(object,"edf")
                  if (is.null(nobs)) nobs <- object@nevent
                  c(-2*logLik(object)/dispersion+k*df+k*df*(df+1)/(nobs-df-1))
              }
          })
setMethod("qAIC", "pstpm2",
          function (object, ..., dispersion = 1, k = 2)  {
              L <- list(...)
              if (length(L)) {
                  L <- c(list(object),L)
                  if (is.null(nobs)) {
                      nobs <- sapply(L,nobs)
                  }
                  if (length(unique(nobs))>1)
                      stop("nobs different: must have identical data for all objects")
                  val <- sapply(L, qAIC, dispersion=dispersion, k=k)
                  df <- sapply(L,attr,"edf")
                  data.frame(qAICc=val,df=df)
              } else {
                  df <- attr(object,"edf")
                  c(-2*logLik(object)/dispersion+k*df)
              }
          })
setMethod("AIC", "pstpm2",
          function (object, ..., k = 2) {
              L <- list(...)
              if (length(L)) {
                  L <- c(list(object),L)
                  if (!all(sapply(L,class)=="pstpm2")) stop("all objects in list must be class pstpm2")
                  val <- sapply(L,AIC,k=k)
                  df <- sapply(L,attr,"edf")
                  data.frame(AIC=val,df=df)
              } else -2 * as.numeric(logLik(object)) + k * attr(object, "edf")
          })
setMethod("BIC", "pstpm2",
          function (object, ..., nobs = NULL) {
              L <- list(...)
              if (length(L)) {
                  L <- c(list(object),L)
                  if (!all(sapply(L,class)=="pstpm2")) stop("all objects in list must be class pstpm2")
                  val <- sapply(L,BIC,nobs=nobs)
                  df <- sapply(L,attr,"edf")
                  data.frame(BIC=val,df=df)
              } else {
                  if (is.null(nobs)) nobs <- object@nevent
                  -2 * as.numeric(logLik(object)) + log(nobs) * attr(object, "edf")
              }
          })
setMethod("eform", signature(object="pstpm2"), 
          function (object, parm, level = 0.95, method = c("Profile"), 
                    name = "exp(beta)") 
          {
              method <- match.arg(method)
              if (missing(parm)) 
                  parm <- TRUE
              estfun <- switch(method, Profile = confint)
              val <- exp(cbind(coef = coef(object), estfun(object, level = level)))
              colnames(val) <- c(name, colnames(val)[-1])
              val[parm, ]
          })
setMethod("show", "pstpm2",
          function(object) {
              object@call.orig <- object@Call
              show(as(object,"mle2"))
              })

simulate.stpm2 <- function(object, nsim=1, seed=NULL,
                           newdata=NULL, lower=1e-6, upper=1e5, start=NULL, ...) {
    if (!is.null(seed)) set.seed(seed)
    if (is.null(newdata)) newdata = as.data.frame(object@data)
    ## assumes nsim replicates per row in newdata
    n = nsim * nrow(newdata)
    if (!is.null(start)) {
        newdatap = newdata
        newdatap[[object@timeVar]] = start # check if this is a sensible size?
        Sentry = predict(object, newdata=newdatap)
        if (length(start)==1)
            lower=rep(start,n)
        else if (length(start)==nrow(newdata))
            lower = rep(start,each=nsim)
        else if (length(start==n))
            lower = start
        else lower = rep(lower,n) # should not get here:(
    } else {
        Sentry = 1
        lower = rep(lower, n)
    }
    newdata = newdata[rep(1:nrow(newdata), each=nsim), , drop=FALSE]
    U <- runif(n)
    objective <- function(time) {
        newdata[[object@timeVar]] <- time
        predict(object, newdata=newdata)/Sentry - U
    }
    vuniroot(objective, lower=rep(lower,length=n), upper=rep(upper,length=n), tol=1e-10, n=n)$root
}
setGeneric("simulate", function(object, nsim=1, seed=NULL, ...) standardGeneric("simulate"))
setMethod("simulate", signature(object="stpm2"),
          function(object, nsim=1, seed=NULL,
                   newdata=NULL, lower=1e-6, upper=1e5, start=NULL, ...)
              simulate.stpm2(object, nsim, seed, newdata, lower,upper,start, ...))
setMethod("simulate", signature(object="pstpm2"), 
          function(object, nsim=1, seed=NULL,
                   newdata=NULL, lower=1e-6, upper=1e5, start=NULL, ...)
              simulate.stpm2(object, nsim, seed, newdata, lower,upper,start, ...))

## Revised from bbmle:
## changed the calculation of the degrees of freedom in the third statement of the .local function
setMethod("anova", signature(object="pstpm2"),
          function (object, ..., width = getOption("width"), 
                    exdent = 10) 
    {
        mlist <- c(list(object), list(...))
        mnames <- sapply(sys.call(sys.parent())[-1], deparse)
        ltab <- as.matrix(do.call("rbind", lapply(mlist, function(x) {
            c(`Tot Df` = x@edf, Deviance = -2 * logLik(x)) # changed to x@edf
        })))
        terms = sapply(mlist, function(obj) {
            if (is.null(obj@formula) || obj@formula == "") {
                mfun <- obj@call$minuslogl
                mfun <- paste("[", if (is.name(mfun)) {
                  as.character(mfun)
                }
                else {
                  "..."
                }, "]", sep = "")
                paste(mfun, ": ", paste(names(obj@coef), collapse = "+"), 
                  sep = "")
            }
            else {
                as.character(obj@formula)
            }
        })
        mterms <- paste("Model ", 1:length(mnames), ": ", mnames, 
            ", ", terms, sep = "")
        mterms <- strwrapx(mterms, width = width, exdent = exdent, 
            wordsplit = "[ \n\t]")
        heading <- paste("Likelihood Ratio Tests", paste(mterms, 
            collapse = "\n"), sep = "\n")
        ltab <- cbind(ltab, Chisq = abs(c(NA, diff(ltab[, "Deviance"]))), 
            Df = abs(c(NA, diff(ltab[, "Tot Df"]))))
        ltab <- cbind(ltab, `Pr(>Chisq)` = c(NA, pchisq(ltab[, 
            "Chisq"][-1], ltab[, "Df"][-1], lower.tail = FALSE)))
        rownames(ltab) <- 1:nrow(ltab)
        attr(ltab, "heading") <- heading
        class(ltab) <- "anova"
        ltab
    })

setMethod("predictnl", "pstpm2",
          function(object,fun,newdata=NULL,link=c("I","log","cloglog","logit"), gd=NULL, ...)
  {
    ## should gd be passed as argument to numDeltaMethod?
    link <- match.arg(link)
    linkf <- eval(parse(text=link))
    if (is.null(newdata) && !is.null(object@data))
      newdata <- object@data
    localf <- function(coef,...)
      {
        object@fullcoef = coef
        linkf(fun(object,...))
      }
    numDeltaMethod(object,localf,newdata=newdata,gd=gd,...)
  })
##
setMethod("predict", "pstpm2",
          function(object,newdata=NULL,
                   type=c("surv","cumhaz","hazard","density","hr","sdiff","hdiff","loghazard","link","meansurv","meansurvdiff","meanhr","odds","or","margsurv","marghaz","marghr","meanhaz","af","fail","margfail","meanmargsurv","rmst","lpmatrix","gradh","gradH","rmstdiff","lpmatrixD"),
                   grid=FALSE,seqLength=300,
                   se.fit=FALSE,link=NULL,exposed=NULL,var=NULL,keep.attributes=FALSE,use.gr=TRUE,level=0.95, n.gauss.quad=100, full=FALSE, ...) {
              type <- match.arg(type)
              predict.stpm2.base(object=object, newdata=newdata, type=type, grid=grid, seqLength=seqLength, se.fit=se.fit,
                                 link=link, exposed=exposed, var=var, keep.attributes=keep.attributes, use.gr=use.gr, level=level, n.gauss.quad=n.gauss.quad, full=full, ...)
              })

setMethod("residuals", "pstpm2",
          function(object, type=c("li","gradli"))
              residuals.stpm2.base(object=object, type=match.arg(type)))

##`%c%` <- function(f,g) function(...) g(f(...)) # function composition
## to do:
## (*) Stata-compatible knots

setMethod("plot", signature(x="pstpm2", y="missing"),
          function(x,y,newdata=NULL,type="surv",
                   xlab=NULL,ylab=NULL,line.col=1,ci.col="grey",lty=par("lty"),
                   add=FALSE,ci=!add,rug=!add,
                   var=NULL,exposed=incrVar(var),times=NULL,...)
              plot.stpm2.base(x=x, y=y, newdata=newdata, type=type, xlab=xlab,
                              ylab=ylab, line.col=line.col, lty=lty, add=add,
                              ci=ci, rug=rug, var=var, exposed=exposed, times=times, ...)
          )
lines.pstpm2 <- function(x,newdata=NULL,type="surv",
                   col=1,ci.col="grey",lty=par("lty"),
                   ci=FALSE,rug=FALSE,
                   var=NULL,exposed=NULL,times=NULL,...)
              plot.stpm2.base(x=x, newdata=newdata, type=type, 
                              line.col=col, ci.col=ci.col, lty=lty, add=TRUE,
                              ci=ci, rug=rug, var=var, exposed=exposed, times=times, ...)
setMethod("lines", signature(x="pstpm2"), lines.pstpm2)

## sandwich variance estimator (from the sandwich package)

## coeftest.stpm2 <- 
## function (x, vcov. = NULL, df = NULL, ...) 
## {
##     est <- coef(x)
##     if (is.null(vcov.)) 
##         se <- vcov(x)
##     else {
##         if (is.function(vcov.)) 
##             se <- vcov.(x)
##         else se <- vcov.
##     }
##     se <- sqrt(diag(se))
##     if (!is.null(names(est)) && !is.null(names(se))) {
##         anames <- names(est)[names(est) %in% names(se)]
##         est <- est[anames]
##         se <- se[anames]
##     }
##     tval <- as.vector(est)/se
##     pval <- 2 * pnorm(abs(tval), lower.tail = FALSE)
##     cnames <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
##     mthd <- "z"
##     rval <- cbind(est, se, tval, pval)
##     colnames(rval) <- cnames
##     class(rval) <- "coeftest"
##     attr(rval, "method") <- paste(mthd, "test of coefficients")
##     return(rval)
## }

## weights.stpm2 <- 
## function (object, ...) 
## {
##     wts <- object@weights
##     if (is.null(wts)) 
##         wts
##     else napredict(object@na.action, wts)
## }

## copy of bbmle:::strwrapx
strwrapx <-
function (x, width = 0.9 * getOption("width"), indent = 0, exdent = 0, 
          prefix = "", simplify = TRUE, parsplit = "\n[ \t\n]*\n", 
          wordsplit = "[ \t\n]") 
{
  if (!is.character(x)) 
    x <- as.character(x)
  indentString <- paste(rep.int(" ", indent), collapse = "")
  exdentString <- paste(rep.int(" ", exdent), collapse = "")
  y <- list()
  plussplit = function(w) {
    lapply(w, function(z) {
      plusloc = which(strsplit(z, "")[[1]] == "+")
      plussplit = apply(cbind(c(1, plusloc + 1), c(plusloc, 
                                                   nchar(z, type = "width"))), 1, function(b) substr(z, 
                                                                                                     b[1], b[2]))
      plussplit
    })
  }
  z <- lapply(strsplit(x, parsplit), function(z) {
    lapply(strsplit(z, wordsplit), function(x) unlist(plussplit(x)))
  })
  for (i in seq_along(z)) {
    yi <- character(0)
    for (j in seq_along(z[[i]])) {
      words <- z[[i]][[j]]
      nc <- nchar(words, type = "w")
      if (any(is.na(nc))) {
        nc0 <- nchar(words)
        nc[is.na(nc)] <- nc0[is.na(nc)]
      }
      if (any(nc == 0)) {
        zLenInd <- which(nc == 0)
        zLenInd <- zLenInd[!(zLenInd %in% (grep("\\.$", 
                                                words) + 1))]
        if (length(zLenInd) > 0) {
          words <- words[-zLenInd]
          nc <- nc[-zLenInd]
        }
      }
      if (length(words) == 0) {
        yi <- c(yi, "", prefix)
        next
      }
      currentIndex <- 0
      lowerBlockIndex <- 1
      upperBlockIndex <- integer(0)
      lens <- cumsum(nc + 1)
      first <- TRUE
      maxLength <- width - nchar(prefix, type = "w") - 
        indent
      while (length(lens) > 0) {
        k <- max(sum(lens <= maxLength), 1)
        if (first) {
          first <- FALSE
          maxLength <- maxLength + indent - exdent
        }
        currentIndex <- currentIndex + k
        if (nc[currentIndex] == 0) 
          upperBlockIndex <- c(upperBlockIndex, currentIndex - 
                                 1)
        else upperBlockIndex <- c(upperBlockIndex, currentIndex)
        if (length(lens) > k) {
          if (nc[currentIndex + 1] == 0) {
            currentIndex <- currentIndex + 1
            k <- k + 1
          }
          lowerBlockIndex <- c(lowerBlockIndex, currentIndex + 
                                 1)
        }
        if (length(lens) > k) 
          lens <- lens[-(1:k)] - lens[k]
        else lens <- NULL
      }
      nBlocks <- length(upperBlockIndex)
      s <- paste(prefix, c(indentString, rep.int(exdentString, 
                                                 nBlocks - 1)), sep = "")
      for (k in (1:nBlocks)) {
        s[k] <- paste(s[k], paste(words[lowerBlockIndex[k]:upperBlockIndex[k]], 
                                  collapse = " "), sep = "")
      }
      s = gsub("\\+ ", "+", s)
      yi <- c(yi, s, prefix)
    }
    y <- if (length(yi)) 
      c(y, list(yi[-length(yi)]))
    else c(y, "")
  }
  if (simplify) 
    y <- unlist(y)
  y
}

## S3 methods
coef.pstpm2 <- coef.stpm2 <- coef
vcov.pstpm2 <- vcov.stpm2 <- vcov

