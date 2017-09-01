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
    q.const <- qr.Q(qr.const, complete=TRUE)[, -(1L:2L), drop = FALSE] # NEW
    basis <- as.matrix((t(qr.qty(qr.const, t(basis))))[, -(1L:nrow(const)), drop = FALSE])
    n.col <- ncol(basis)
    if (nas) {
        nmat <- matrix(NA, length(nax), n.col)
        nmat[!nax, ] <- basis
        basis <- nmat
    }
    dimnames(basis) <- list(nx, 1L:n.col)
    if (centre) {
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
grad1 <- function(func,x,...)
  {
    h <- .Machine$double.eps^(1/3)*ifelse(abs(x)>1,abs(x),1)
    temp <- x+h
    h.hi <- temp-x
    temp <- x-h
    h.lo <- x-temp
    twoeps <- h.hi+h.lo
    ny <- length(func(x,...))
    if (ny==0L) stop("Length of function equals 0")
    (func(x+h, ...) - func(x-h, ...))/twoeps
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
numDeltaMethod <- function(object,fun,gd=NULL,...) {
  coef <- coef(object)
  est <- fun(coef,...)
  Sigma <- vcov(object)
  if (is.null(gd))
      gd <- grad(fun,coef,...)
  ## se.est <- as.vector(sqrt(diag(t(gd) %*% Sigma %*% gd)))
  se.est <- as.vector(sqrt(colSums(gd* (Sigma %*% gd))))
  data.frame(Estimate = est, SE = se.est)
}
"coef<-" <- function (x, value) 
  UseMethod("coef<-")
predictnl <- function (object, ...) 
  UseMethod("predictnl")
"coef<-.default" <- function(x,value) {
    x$coefficients <- value
    x
}
predictnl.default <- function(object,fun,newdata=NULL,gd=NULL,...)
  {
    ## link=c(I,log,sqrt),invlink=NULL
    ## link <- match.arg(link)
    ## if (is.null(invlink))
    ##       invlink <- switch(deparse(substitute(link)),I=I,log=exp,sqrt=function(x) x^2)
    if (is.null(newdata) && !is.null(object$data))
      newdata <- object$data
    localf <- function(coef,...)
      {
        if ("coefficients" %in% names(object)) {
            object$coefficients <- coef
        } else if ("coef" %in% names(object)) {
            object$coef <- coef
        } else coef(object) <- coef
        fun(object,...)
      }
    numDeltaMethod(object,localf,newdata=newdata,gd=gd,...)
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
## setMethod("predictnl", "mle", function(object,fun,gd=NULL,...)
##   {
##     localf <- function(coef,...)
##       {
##         object@fullcoef = coef # changed from predictnl.default()
##         fun(object,...)
##       }
##     numDeltaMethod(object,localf,gd=gd,...)
##   })
predict.formula <- function(formula,data,newdata,na.action,type="model.matrix") 
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
stpm2 <- function(formula, data, smooth.formula = NULL, smooth.args = NULL,
                     df = 3, cure = FALSE, logH.args = NULL, logH.formula = NULL,
                     tvc = NULL, tvc.formula = NULL,
                     control = list(parscale = 1, maxit = 300), init = NULL,
                     coxph.strata = NULL, weights = NULL, robust = FALSE, baseoff = FALSE, 
                  bhazard = NULL, timeVar = "", time0Var = "", use.gr = TRUE,
                  optimiser=c("BFGS","NelderMead"),
                     reltol=1.0e-8, trace = 0,
                     link.type=c("PH","PO","probit","AH","AO"), theta.AO=0, 
                  frailty = !is.null(cluster) & !robust, cluster = NULL, logtheta=-6, nodes=9, RandDist=c("Gamma","LogN"), recurrent = FALSE,
                  adaptive = TRUE, maxkappa = 1e3, Z = ~1,
                     contrasts = NULL, subset = NULL, robust_initial=FALSE, ...) {
    link.type <- match.arg(link.type)
    link <- switch(link.type,PH=link.PH,PO=link.PO,probit=link.probit,AH=link.AH,AO=link.AO(theta.AO))
    RandDist <- match.arg(RandDist)
    optimiser <- match.arg(optimiser)
    use.gr <- TRUE # old code
    ## logH.formula and logH.args are deprecated
    if (!is.null(smooth.formula) && is.null(logH.formula))
        logH.formula <- smooth.formula
    if (!is.null(smooth.args) && is.null(logH.args))
        logH.args <- smooth.args
    ## parse the event expression
    eventInstance <- eval(lhs(formula),envir=data)
    stopifnot(length(lhs(formula))>=2)
    eventExpr <- lhs(formula)[[length(lhs(formula))]]
    delayed <- length(lhs(formula))>=4 # indicator for multiple times (cf. strictly delayed)
    surv.type <- attr(eventInstance,"type")
    if (surv.type %in% c("interval2","left","mstate"))
        stop("stpm2 not implemented for Surv type ",surv.type,".")
    counting <- attr(eventInstance,"type") == "counting"
    interval <- attr(eventInstance,"type") == "interval"
    timeExpr <- lhs(formula)[[if (delayed) 3 else 2]] # expression
    if (timeVar == "")
        timeVar <- all.vars(timeExpr)
    ## set up the formulae
    if (is.null(logH.formula) && is.null(logH.args)) {
        logH.args$df <- df
        if (cure) logH.args$cure <- cure
    }
    if (is.null(logH.formula))
      logH.formula <- as.formula(call("~",as.call(c(quote(nsx),call("log",timeExpr),
                                                    vector2call(logH.args)))))
    if (is.null(tvc.formula) && !is.null(tvc)) {
      tvc.formulas <-
        lapply(names(tvc), function(name)
               call(":",
                    as.name(name),
                    as.call(c(quote(nsx),
                              call("log",timeExpr),
                              vector2call(if (cure) list(cure=cure,df=tvc[[name]]) else list(df=tvc[[name]])
                                          )))))
      if (length(tvc.formulas)>1)
        tvc.formulas <- list(Reduce(`%call+%`, tvc.formulas))
      tvc.formula <- as.formula(call("~",tvc.formulas[[1]]))
    }
    if (!is.null(tvc.formula)) {
      rhs(logH.formula) <- rhs(logH.formula) %call+% rhs(tvc.formula)
    }
    if (baseoff)
      rhs(logH.formula) <- rhs(tvc.formula)
    full.formula <- formula
    rhs(full.formula) <- rhs(formula) %call+% rhs(logH.formula)
    Z.formula <- Z
    ##
    ## set up the data
    ## ensure that data is a data frame
    ## data <- get_all_vars(full.formula, data) # but this loses the other design information
    ## restrict to non-missing data (assumes na.action=na.omit)
    .include <- apply(model.matrix(formula, data, na.action = na.pass), 1, function(row) !any(is.na(row))) &
        !is.na(eval(eventExpr,data)) & !is.na(eval(timeExpr,data))
    data <- data[.include, , drop=FALSE]
    ##
    ## parse the function call
    Call <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "contrasts", "weights"),
               names(mf), 0L)
    mf <- mf[c(1L, m)]
    ##
    ## get variables
    time <- eval(timeExpr, data, parent.frame())
    if (any(time>0 & time<1e-4))
        warning("Some event times < 1e-4: consider transforming time to avoid problems with finite differences")
    time0Expr <- NULL # initialise
    if (delayed) {
      time0Expr <- lhs(formula)[[2]]
      if (time0Var == "")
        time0Var <- all.vars(time0Expr)
      time0 <- eval(time0Expr, data, parent.frame())
      if (any(time0>0 & time0<1e-4))
          warning("Some entry times < 1e-4: consider transforming time to avoid problems with finite differences")
    }
    event <- eval(eventExpr,data)
    ## if all the events are the same, we assume that they are all events, else events are those greater than min(event)
    event <- if (length(unique(event))==1) rep(TRUE, length(event)) else event <- event > min(event)
    ## setup for initial values
    if (!interval) {
        ## Cox regression
        coxph.call <- mf
        coxph.call[[1L]] <- as.name("coxph")
        coxph.strata <- substitute(coxph.strata)
        if (!is.null(coxph.strata)) {
            coxph.formula <- formula
            rhs(coxph.formula) <- rhs(formula) %call+% call("strata",coxph.strata)
            coxph.call$formula <- coxph.formula
        }
        coxph.call$model <- TRUE
        coxph.obj <- eval(coxph.call, envir=parent.frame())
        y <- model.extract(model.frame(coxph.obj),"response")
        data$logHhat <- if (is.null(bhazard)) {
                            pmax(-18,link$link(Shat(coxph.obj)))
                        } else  pmax(-18,link$link(Shat(coxph.obj)/exp(-bhazard*time)))
    }
    if (interval) {
        ## survref regression
        survreg.call <- mf
        survreg.call[[1L]] <- as.name("survreg")
        survreg.obj <- eval(survreg.call, envir=parent.frame())
        weibullShape <- 1/survreg.obj$scale
        weibullScale <- predict(survreg.obj)
        y <- model.extract(model.frame(survreg.obj),"response")
        data$logHhat <- pmax(-18,link$link(pweibull(time,weibullShape,weibullScale,lower.tail=FALSE)))
    }
    ##
    ## initial values and object for lpmatrix predictions
    lm.call <- mf
    lm.call[[1L]] <- as.name("lm")
    lm.formula <- full.formula
    lhs(lm.formula) <- quote(logHhat) # new response
    lm.call$formula <- lm.formula
    dataEvents <- data[event,]
    if (interval)
        dataEvents <- data
    lm.call$data <- quote(dataEvents) # events only
    lm.obj <- eval(lm.call)
    if (is.null(init)) {
      init <- coef(lm.obj)
    }
    ##
    ## set up mf and wt
    mt <- terms(lm.obj)
    mf <- model.frame(lm.obj)
    ## wt <- model.weights(lm.obj$model)
    wt <- if (is.null(substitute(weights))) rep(1,nrow(data)) else eval(substitute(weights),data,parent.frame())
    ##
    ## XD matrix
    lpfunc <- function(delta,fit,dataset,var) {
      dataset[[var]] <- dataset[[var]]+delta
      lpmatrix.lm(fit,dataset)
    }
    ##
    bhazard <- substitute(bhazard)
    bhazard <- if (is.null(bhazard)) rep(0,nrow(data)) else eval(bhazard,data,parent.frame())
    excess <- !all(bhazard==0)
    ## initialise values specific to either delayed entry or interval-censored
    ind0 <- FALSE
    map0 <- 0L
    which0 <- 0
    wt0 <- 0
    ttype <- 0
    transX <- function(X, data) X
    transXD <- function(XD) XD
    if (!interval) { # surv.type %in% c("right","counting")
        X <- lpmatrix.lm(lm.obj,data)
        if (link.type=="AH") {
            datat0 <- data
            datat0[[timeVar]] <- 0
            index0 <- which.dim(X - lpmatrix.lm(lm.obj,datat0))
            transX <- function(X, data) {
                datat0 <- data
                datat0[[timeVar]] <- 0
                Xt0 <- lpmatrix.lm(lm.obj,datat0)
                (X - Xt0)[, index0, drop=FALSE]
            }
            transXD <- function(XD) XD[, index0, drop=FALSE]
            init <- init[index0]
        }
        X <- transX(X,data)
        XD <- grad(lpfunc,0,lm.obj,data,timeVar)
        XD <- transXD(matrix(XD,nrow=nrow(X)))
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
            X0 <- transX(lpmatrix.lm(lm.obj, data0), data0)
            wt0 <- wt[ind0]
            rm(data0)
        }
    } else { ## interval-censored
        ## ttime <- eventInstance[,1]
        ## ttime2 <- eventInstance[,2]
        ttype <- eventInstance[,3]
        X1 <- transX(lpmatrix.lm(lm.obj,data),data)
        data0 <- data
        data0[[timeVar]] <- data0[[time0Var]]
        X <- transX(lpmatrix.lm(lm.obj, data0), data0)
        XD <- grad(lpfunc,0,lm.obj,data0,timeVar)
        XD <- transXD(matrix(XD,nrow=nrow(X)))
        X0 <- matrix(0,nrow(X),ncol(X))
        rm(data0)
    } 
    if (frailty) {
        Z <- model.matrix(Z, data)
        if (ncol(Z)>2) stop("Current implementation only allows for one or two random effects")
        if (ncol(Z)==2) {
            init <- c(init,logtheta1=logtheta,corrtrans=0,logtheta2=logtheta)
        } else init <- c(init, logtheta=logtheta)
    }
    if (!is.null(control) && "parscale" %in% names(control)) {
      if (length(control$parscale)==1)
        control$parscale <- rep(control$parscale,length(init))
      if (is.null(names(control$parscale)))
        names(control$parscale) <- names(init)
    }
    parscale <- rep(if (is.null(control$parscale)) 1 else control$parscale,length=length(init))
    names(parscale) <- names(init)
    args <- list(init=init,X=X,XD=XD,bhazard=bhazard,wt=wt,event=ifelse(event,1,0),time=time,
                 delayed=delayed, interval=interval, X0=X0, wt0=wt0, X1=X1, parscale=parscale, reltol=reltol,
                 kappa=1, trace = trace, oldcluster=cluster, frailty=frailty, cluster=if(!is.null(cluster)) as.vector(unclass(factor(cluster))) else NULL, map0 = map0 - 1L, ind0 = ind0, which0 = which0 - 1L, link=link.type, ttype=ttype,
                 RandDist=RandDist, optimiser=optimiser,
                 type=if (frailty && RandDist=="Gamma") "stpm2_gamma_frailty" else if (frailty && RandDist=="LogN") "stpm2_normal_frailty" else "stpm2", recurrent = recurrent, return_type="optim", transX=transX, transXD=transXD, maxkappa=maxkappa, Z=Z, Z.formula = Z.formula, thetaAO = theta.AO, excess=excess, data=data,
                 robust_initial = robust_initial)
    if (frailty) {
        rule <- fastGHQuad::gaussHermiteData(nodes)
        args$gauss_x <- rule$x
        args$gauss_w <- rule$w
        args$adaptive <- adaptive
        if (ncol(args$Z)>1) {
            use.gr <- FALSE
            args$type <- "stpm2_normal_frailty_2d"
            args$Z <- as.matrix(args$Z)
            ## args$adaptive <- FALSE # adaptive not defined
            ## args$optimiser <- "NelderMead" # gradients not defined
        } else args$Z <- as.vector(args$Z)
    }
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
    parnames(negll) <- parnames(gradnegll) <- names(init)
    ## MLE
    fit <- .Call("model_output", args, PACKAGE="rstpm2")
    args$init <- coef <- as.vector(fit$coef)
    args$kappa.final <- fit$kappa
    hessian <- fit$hessian
    names(coef) <- rownames(hessian) <- colnames(hessian) <- names(init)
    mle2 <- if (use.gr) mle2(negll, coef, vecpar=TRUE, control=control, gr=gradnegll, ..., eval.only=TRUE) else mle2(negll, coef, vecpar=TRUE, control=control, ..., eval.only=TRUE)
    mle2@vcov <- if (!inherits(vcov <- try(solve(hessian)), "try-error")) vcov else matrix(NA,length(coef), length(coef))
    mle2@details$convergence <- fit$fail # fit$itrmcd
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
               data = data,
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
    return(out)
  }
## summary.mle is not exported from bbmle
.__C__summary.mle2 <- bbmle:::.__C__summary.mle2 # hack suggested from http://stackoverflow.com/questions/28871632/how-to-resolve-warning-messages-metadata-object-not-found-spatiallinesnull-cla
setClass("summary.stpm2", representation(frailty="logical",theta="list",wald="matrix",args="list"), contains="summary.mle2")
## setAs("summary.stpm2", "summary.mle2",
##       function(from,to) new("summary.mle2", call=from@call, coef=from@call, m2logL=from@m2logL))
## setMethod("show", "stpm2", function(object) show(as(object,"mle2")))
corrtrans <- function(x) (1-exp(-x)) / (1+exp(-x))
setMethod("summary", "stpm2",
          function(object) {
              newobj <- as(summary(as(object,"mle2")),"summary.stpm2")
              newobj@args <- object@args
              newobj@frailty <- object@frailty
              if (object@frailty && !is.matrix(object@args$Z)) {
                  coef <- coef(newobj)
                  theta <- exp(coef[nrow(coef),1])
                  se.logtheta <- coef[nrow(coef),2]
                  se.theta <- theta*se.logtheta
                  test.statistic <- (1/se.logtheta)^2
                  p.value <- pchisq(test.statistic,df=1,lower.tail=FALSE)/2
                  newobj@theta <- list(theta=theta, se.theta=se.theta, p.value=p.value)
              } else if (object@frailty && is.matrix(object@args$Z) && ncol(object@args$Z)==2) {
                  ## browser()
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
              if (object@frailty) {
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
    invlinkf <- switch(link,I=I,log=exp,cloglog=cexpexp,logit=expit)
    linkf <- eval(parse(text=link))
    if (is.null(newdata) && !is.null(object@data))
      newdata <- object@data
    localf <- function(coef,...)
      {
        object@fullcoef = coef
        linkf(fun(object,...))
      }
    dm <- numDeltaMethod(object,localf,newdata=newdata,gd=gd,...)
    out <- invlinkf(data.frame(Estimate=dm$Estimate,
                               lower=dm$Estimate-1.96*dm$SE,
                               upper=dm$Estimate+1.96*dm$SE))
    ## cloglog switches the bounds
    if (link=="cloglog") 
      out <- data.frame(Estimate=out$Estimate,lower=out$upper,upper=out$lower)
    return(out)
  })
##
setMethod("predictnl", "aft",
          function(object,fun,newdata=NULL,link=c("I","log","cloglog","logit"), gd=NULL, ...)
  {
    link <- match.arg(link)
    invlinkf <- switch(link,I=I,log=exp,cloglog=cexpexp,logit=expit)
    linkf <- eval(parse(text=link))
    if (is.null(newdata) && !is.null(object@args$data))
      newdata <- object@args$data
    localf <- function(coef,...)
      {
          object@args$init <- object@fullcoef <- coef
        linkf(fun(object,...))
      }
    dm <- numDeltaMethod(object,localf,newdata=newdata,gd=gd,...)
    out <- invlinkf(data.frame(Estimate=dm$Estimate,
                               lower=dm$Estimate-1.96*dm$SE,
                               upper=dm$Estimate+1.96*dm$SE))
    ## cloglog switches the bounds
    if (link=="cloglog") 
      out <- data.frame(Estimate=out$Estimate,lower=out$upper,upper=out$lower)
    return(out)
  })

residuals.stpm2.base <- function(object, type=c("li","gradli")) {
    type <- match.arg(type)
    args <- object@args
    if (type=="li") {
        localargs <- args
        localargs$return_type <- "li"
        return(as.vector(.Call("model_output", localargs, PACKAGE="rstpm2")))
    }
    if (type=="gradli") {
        localargs <- args
        localargs$return_type <- "gradli"
        return(.Call("model_output", localargs, PACKAGE="rstpm2"))
    }
}    
setMethod("residuals", "stpm2",
          function(object, type=c("li","gradli"))
              residuals.stpm2.base(object=object, type=type))   

predict.stpm2.base <- 
          function(object,newdata=NULL,
                   type=c("surv","cumhaz","hazard","density","hr","sdiff","hdiff","loghazard","link","meansurv","meansurvdiff","odds","or","margsurv","marghaz","marghr","meanhaz","af","fail","margfail","meanmargsurv","uncured"),
                   grid=FALSE,seqLength=300,
                   se.fit=FALSE,link=NULL,exposed=NULL,var=NULL,keep.attributes=TRUE,use.gr=TRUE,...)
{
    type <- match.arg(type)
    args <- object@args
    if (type %in% c("fail","margfail")) {
        out <- 1-predict.stpm2.base(object,newdata=newdata,type=switch(type, fail="surv", margfail="margsurv"),grid,seqLength,se.fit,link,exposed,var,keep.attributes,use.gr,...)
        if (se.fit) {temp <- out$lower; out$lower <- out$upper; out$upper <- temp}
        return(out)
    }
      if (is.null(exposed) && is.null(var) & type %in% c("hr","sdiff","hdiff","meansurvdiff","or","marghr","af","uncured"))
          stop('Either exposed or var required for type in ("hr","sdiff","hdiff","meansurvdiff","or","marghr","af","uncured")')
      if (type %in% c('margsurv','marghaz','marghr','margfail','meanmargsurv') && !object@args$frailty)
          stop("Marginal prediction only for frailty models")
    ## exposed is a function that takes newdata and returns the revised newdata
    ## var is a string for a variable that defines a unit change in exposure
    if (is.null(newdata) && type %in% c("hr","sdiff","hdiff","meansurvdiff","or","marghr","uncured"))
        stop("Prediction using type in ('hr','sdiff','hdiff','meansurvdiff','or','marghr','uncured') requires newdata to be specified.")
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
    lpfunc <- if (inherits(object,"pstpm2")) 
        function(x,...) {
            newdata2 <- newdata
            newdata2[[object@timeVar]] <- x
            predict(object@gam,newdata2,type="lpmatrix")
        } else  function(delta,fit,data,var) {
            data[[var]] <- data[[var]]+delta
            lpmatrix.lm(fit,data)
        }
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
    if (calcX)  {
      if (inherits(object, "stpm2")) {
          X <- object@args$transX(lpmatrix.lm(object@lm, newdata), newdata)
          XD <- grad(lpfunc,0,object@lm,newdata,object@timeVar)
          XD <- object@args$transXD(matrix(XD,nrow=nrow(X)))
      }
      if (inherits(object, "pstpm2")) {
           X <- object@args$transX(predict(object@gam, newdata, type="lpmatrix"), newdata)      
           XD <- object@args$transXD(grad1(lpfunc,newdata[[object@timeVar]]))
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
    ## if (object@delayed && !object@interval) {
    ##   newdata0 <- newdata
    ##   newdata0[[object@timeVar]] <- newdata[[object@time0Var]]
    ##   X0 <- lpmatrix.lm(object@lm, newdata0)
    ##   ## XD0 <- grad(lpfunc,0,object@lm,newdata,object@timeVar)
    ##   ## XD0 <- matrix(XD0,nrow=nrow(X0))
    ## }
    if (type %in% c("hr","sdiff","hdiff","meansurvdiff","or","marghr","af","uncured")) {
        newdata2 <- exposed(newdata)
      if (inherits(object, "stpm2")) {
          X2 <- object@args$transX(lpmatrix.lm(object@lm, newdata2), newdata2)
          XD2 <- grad(lpfunc,0,object@lm,newdata2,object@timeVar)
          XD2 <- object@args$transXD(matrix(XD2,nrow=nrow(X)))
      }
      if (inherits(object, "pstpm2")) {
           X2 <- object@args$transX(predict(object@gam, newdata2, type="lpmatrix"), newdata2)      
           XD2 <- object@args$transXD(grad1(lpfunc,newdata2[[object@timeVar]]))
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
        link <- object@link
        beta <- coef(object)
        npar <- length(beta)
        logtheta <- beta[npar]
        theta <- exp(beta[npar])
        beta <- beta[-npar]
        Hessian <- solve(vcov(object))
        eta <- as.vector(X %*% beta)
        eta2 <- as.vector(X2 %*% beta)
        S <- link$ilink(eta)
        S2 <- link$ilink(eta2)
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
                               rbind(colSums(margS[index]*(-link$gradH(eta[index],list(X=X[index,,drop=FALSE]))/(1+theta*H[index])))/n.cluster,
                                     colSums(margS2[index]*(-link$gradH(eta2[index],list(X=X2[index,,drop=FALSE]))/(1+theta*H2[index])))/n.cluster),
                               c(sum(dmarg.dlogtheta(logtheta,H[index]))/n.cluster,
                                 sum(dmarg.dlogtheta(logtheta,H2[index]))/n.cluster))
            par.hessian <- cbind(matrix(0, nrow = npar, ncol = 2), -Hessian / n.cluster)
            bread <- rbind(S.hessian, par.hessian)
            ibread <- solve(bread)
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
        ## browser()
        times <- newdata[[object@timeVar]]
        utimes <- sort(unique(times))
        n <- nrow(X)/length(utimes)
        n.cluster <- length(unique(args$cluster))
        link <- object@link
        beta <- coef(object)
        npar <- length(beta)
        logtheta <- beta[npar]
        theta <- exp(beta[npar])
        beta <- beta[-npar]
        Hessian <- solve(vcov(object))
        eta <- as.vector(X %*% beta)
        S <- link$ilink(eta)
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
            ## browser()
            index <- which(times==utimes[i])
            newobj <- object
            newobj@args$X <- X[index,,drop=FALSE]
            newobj@args$XD <- XD[index,,drop=FALSE]
            ## ## Attempt to use the sandwich estimator for the covariance matrix with the delta method (-> inflated SEs)
            ## res <- residuals(newobj,type="gradli")
            ## meat <- stats::var(res, na.rm=TRUE) # crossprod(res)/n.cluster
            ## colnames(meat) <- rownames(meat) <- c(names(beta),"logtheta")
            ## bread <- -Hessian / n.cluster
            ## ibread <- solve(bread)
            ## sandwich <- ibread %*% meat %*% t(ibread) / n.cluster
            ## g <- c(colSums(margS[index]*(-link$gradH(eta[index],list(X=X[index,,drop=FALSE]))/(1+theta*H[index])))/n.cluster,
            ##                    sum(dmarg.dlogtheta(logtheta,H[index]))/n.cluster)
            ## se.fit[i] <- sqrt(sum(matrix(g,nrow=1)%*%(sandwich %*% g)))
            gradli <- residuals(newobj,type="gradli")
            res <- tapply(margS[index]-mean(margS[index]),args$cluster,sum)
            res <- cbind(res, gradli)
            meat <- stats::var(res, na.rm=TRUE)
            ## meat <- crossprod(res)/n.cluster
            colnames(meat) <- rownames(meat) <- c("S", names(beta),"logtheta")
            S.hessian <- c(-n/n.cluster,
                               colSums(margS[index]*(-link$gradH(eta[index],list(X=X[index,,drop=FALSE]))/(1+theta*H[index])))/n.cluster,
                               sum(dmarg.dlogtheta(logtheta,H[index]))/n.cluster)
            par.hessian <- cbind(matrix(0, nrow = npar, ncol = 1), -Hessian / n.cluster)
            bread <- rbind(S.hessian, par.hessian)
            ibread <- solve(bread)
            sandwich <- (ibread %*% meat %*% t(ibread) / n.cluster)[1, 1]
            se.fit[i] <- sqrt(sandwich)
        }
        pred <- data.frame(Estimate=fit, lower=fit-1.96*se.fit, upper=fit+1.96*se.fit) 
        if (keep.attributes)
            attr(pred,"newdata") <- newdata
        return(pred)
    }
    local <-  function (object, newdata=NULL, type="surv", exposed)
    {
        beta <- coef(object)
        tt <- object@terms
        link <- object@link
        if (object@frailty) {
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
        if (!object@args$excess && any(h<0)) warning(sprintf("Predicted hazards less than zero (n=%i).",sum(h<0)))
        H = link$H(eta)
        Sigma = vcov(object)
        if (type=="link") {
          return(eta)
        }
        if (type=="cumhaz") {
            ## if (object@delayed) {
            ##     eta0 <- as.vector(X0 %*% beta)
            ##     etaD0 <- as.vector(XD0 %*% beta)
            ##     H0 <- link$H(eta0)
            ##     return(H - H0)
            ## }
            ## else 
                return(H)
        }
        if (type=="density")
            return (S*h)
        if (type=="surv") {
          return(S)
        }
        if (type=="fail") {
          return(1-S)
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
    }
    pred <- if (!se.fit) {
      local(object,newdata,type=type,exposed=exposed,
            ...)
    }
    else {
      gd <- NULL
      beta <- coef(object)
      if (is.null(link))
        link <- switch(type,surv="cloglog",cumhaz="log",hazard="log",hr="log",sdiff="I",
                       hdiff="I",loghazard="I",link="I",odds="log",or="log",margsurv="cloglog",marghaz="log",marghr="log",meansurv="I",meanhaz="I",af="I",uncured="cloglog")
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
          if (type=="hazard" && link %in% c("I","log")) {
              ## Case: frailty model (assumes baseline hazard for frailty=1)
              betastar <- if(args$frailty) beta[-length(beta)] else beta
              gd <- switch(link,
                           I=t(object@link$gradh(X %*% betastar, XD %*% betastar, list(X=X, XD=XD))),
                           log=t(object@link$gradh(X %*% betastar, XD %*% betastar, list(X=X, XD=XD))/
                               object@link$h(X %*% betastar, XD %*% betastar)))
          }
          if (type=="meansurv" && !object@frailty) {
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
      predictnl(object,local,link=link,newdata=newdata,type=type,gd=if (use.gr) gd else NULL,
                exposed=exposed,...) 
    }
    if (keep.attributes)
        attr(pred,"newdata") <- newdata
    return(pred)
  }

setMethod("predict", "stpm2",
          function(object,newdata=NULL,
                   type=c("surv","cumhaz","hazard","density","hr","sdiff","hdiff","loghazard","link","meansurv","meansurvdiff","odds","or","margsurv","marghaz","marghr","meanhaz","af","fail","margfail","meanmargsurv","uncured"),
                   grid=FALSE,seqLength=300,
                   se.fit=FALSE,link=NULL,exposed=incrVar(var),var=NULL,keep.attributes=TRUE,use.gr=TRUE,...)
              predict.stpm2.base(object=object, newdata=newdata, type=type, grid=grid, seqLength=seqLength, se.fit=se.fit,
                                 link=link, exposed=exposed, var=var, keep.attributes=keep.attributes, use.gr=use.gr, ...))

##`%c%` <- function(f,g) function(...) g(f(...)) # function composition
plot.meansurv <- function(x, y=NULL, times=NULL, newdata=NULL, type="meansurv", exposed=NULL, add=FALSE, ci=!add, rug=!add, recent=FALSE,
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
        pred <- predict(x, newdata=newdata, type=type, se.fit=ci, exposed=exposed) # requires recent version
        if (type=="meansurv")
            pred <- if (ci) rbind(c(Estimate=1,lower=1,upper=1),pred) else c(1,pred)
    } else {
        pred <- lapply(times, 
                       function(time) {
                           newdata[[x@timeVar]] <- newdata[[x@timeVar]]*0+time
                           predict(x, newdata=newdata, type=type, se.fit=ci, grid=FALSE, exposed=exposed)
                       })
        pred <- do.call("rbind", pred)
        if (type=="meansurv")  {
            pred <- if (ci) rbind(c(Estimate=1,lower=1,upper=1),pred) else c(1,unlist(pred))
            times <- c(0,times)
            }
        }
    if (is.null(xlab)) xlab <- deparse(x@timeExpr)
    if (is.null(ylab))
        ylab <- switch(type,
                       meansurv="Mean survival",
                       af="Attributable fraction",
                       meansurvdiff="Difference in mean survival")
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
                   xlab=NULL,ylab=NULL,line.col=1,ci.col="grey",lty=par("lty"),
                   add=FALSE,ci=!add,rug=!add,
                   var=NULL,exposed=incrVar(var),times=NULL,...) {
              if (type %in% c("meansurv","meansurvdiff","af")) {
                  return(plot.meansurv(x,times=times,newdata=newdata,type=type,xlab=xlab,ylab=ylab,line.col=line.col,ci.col=ci.col,
                                       lty=lty,add=add,ci=ci,rug=rug, exposed=exposed, ...))
              }
              if (is.null(newdata)) stop("newdata argument needs to be specified")
              y <- predict(x,newdata,type=switch(type,fail="surv",margfail="margsurv",type),var=var,exposed=exposed,
                           grid=!(x@timeVar %in% names(newdata)), se.fit=ci)
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
                                 meansurvdiff="Difference in mean survival",odds="Odds",or="Odds ratio",
                                 margsurv="Marginal survival",marghaz="Marginal hazard",marghr="Marginal hazard ratio", haz="Hazard",fail="Failure",
                                 meanhaz="Mean hazard",margfail="Marginal failure",af="Attributable fraction",meanmargsurv="Mean marginal survival",
                                 uncured="Uncured distribution")
              xx <- attr(y,"newdata")
              xx <- eval(x@timeExpr,xx) # xx[,ncol(xx)]
              if (!add) matplot(xx, y, type="n", xlab=xlab, ylab=ylab, ...)
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
                   var=NULL,exposed=incrVar(var),times=NULL,...)
              plot.stpm2.base(x=x, y=y, newdata=newdata, type=type, xlab=xlab,
                              ylab=ylab, line.col=line.col, lty=lty, add=add,
                              ci=ci, rug=rug, var=var, exposed=exposed, times=times, ...)
          )
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
	                          fullformula="formula",
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
pstpm2 <- function(formula, data, smooth.formula = NULL, smooth.args = NULL,
                   logH.args = NULL, 
                   tvc = NULL, 
                   control = list(parscale = 1, maxit = 300), init = NULL,
                   coxph.strata = NULL, coxph.formula = NULL,
                   weights = NULL, robust = FALSE, 
                   bhazard = NULL, timeVar = "", time0Var = "",
                   sp=NULL, use.gr = TRUE, 
                   criterion=c("GCV","BIC"), penalty = c("logH","h"), smoother.parameters = NULL,
                   alpha=if (is.null(sp)) switch(criterion,GCV=1,BIC=1) else 1, sp.init=1, trace = 0,
                   link.type=c("PH","PO","probit","AH","AO"), theta.AO=0,
                  optimiser=c("BFGS","NelderMead","Nlm"), recurrent = FALSE,
                  frailty=!is.null(cluster) & !robust, cluster = NULL, logtheta=-6, nodes=9,RandDist=c("Gamma","LogN"),
                  adaptive=TRUE, maxkappa = 1e3, Z = ~1,
                   reltol = list(search = 1.0e-10, final = 1.0e-10, outer=1.0e-5),outer_optim=1,
                   contrasts = NULL, subset = NULL, robust_initial = FALSE, ...) {
    link.type <- match.arg(link.type)
    link <- switch(link.type,PH=link.PH,PO=link.PO,probit=link.probit,AH=link.AH,AO=link.AO(theta.AO))
    RandDist <- match.arg(RandDist)
    optimiser <- match.arg(optimiser)
    ## logH.args is deprecated
    if (!is.null(smooth.args) && is.null(logH.args))
        logH.args <- smooth.args
    ## set up the data
    ## ensure that data is a data frame
    temp.formula <- formula
    if (!is.null(smooth.formula)) rhs(temp.formula) <-rhs(temp.formula) %call+% rhs(smooth.formula)
    raw.data <- data
    ## data <- get_all_vars(temp.formula, raw.data)
    criterion <- match.arg(criterion)
    penalty <- match.arg(penalty)
    ##
    ## parse the function call
    Call <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "contrasts", "weights"),
               names(mf), 0L)
    mf <- mf[c(1L, m)]
    ##
    ## parse the event expression
    eventExpr <- lhs(formula)[[length(lhs(formula))]]
    eventInstance <- eval(lhs(formula),envir=data)
    stopifnot(length(lhs(formula))>=2)
    delayed <- length(lhs(formula))>=4
    surv.type <- attr(eventInstance,"type")
    if (surv.type %in% c("interval2","left","mstate"))
        stop("stpm2 not implemented for Surv type ",surv.type,".")
    interval <- attr(eventInstance,"type") == "interval"
    timeExpr <- lhs(formula)[[if (delayed) 3 else 2]] # expression
    if (timeVar == "")
      timeVar <- all.vars(timeExpr)
    ## restrict to non-missing data (assumes na.action=na.omit)
    .include <- apply(model.matrix(formula, data, na.action = na.pass), 1, function(row) !any(is.na(row))) &
        !is.na(eval(eventExpr,data)) & !is.na(eval(timeExpr,data))
    data <- data[.include, , drop=FALSE] ### REPLACEMENT ###
    ## we can now evaluate over data
    time <- eval(timeExpr, data, parent.frame())
    time0Expr <- NULL # initialise
    if (delayed) {
      time0Expr <- lhs(formula)[[2]]
      if (time0Var == "")
        time0Var <- all.vars(time0Expr)
      time0 <- eval(time0Expr, data, parent.frame())
    }
    event <- eval(eventExpr,data,parent.frame())
    event <- event > min(event)
    nevent <- sum(event)
    ##
    ## set up the formulae
    if (is.null(smooth.formula) && is.null(logH.args)) {
      logH.args$k <- -1
    }
    if (is.null(smooth.formula))
      smooth.formula <- as.formula(call("~",as.call(c(quote(s),call("log",timeExpr),
                                                    vector2call(logH.args)))))
    if (!is.null(tvc)) {
      tvc.formulas <-
        lapply(names(tvc), function(name)
               call(":",
                    as.name(name),
                    as.call(c(quote(s),
                              call("log",timeExpr),
                              vector2call(list(k=tvc[[name]]))))))
      if (length(tvc.formulas)>1)
        tvc.formulas <- list(Reduce(`%call+%`, tvc.formulas))
      tvc.formula <- as.formula(call("~",tvc.formulas[[1]]))
      rhs(smooth.formula) <- rhs(smooth.formula) %call+% rhs(tvc.formula)
    }
    full.formula <- formula
    if(link.type=="AH"){
      rhs(full.formula) <- rhs(smooth.formula)
    }
    else{
      rhs(full.formula) <- rhs(formula) %call+% rhs(smooth.formula)
    }
      ## 
	  left <- deparse(substitute(formula))
	  tf <- terms.formula(smooth.formula, specials = c("s", "te"))
		terms <- attr(tf, "term.labels")
		right <- paste0(terms, collapse = "+")
    fullformula <- as.formula(paste0(left, "+", right), env = parent.frame())
    if (!interval) {
        ## Cox regression
        coxph.call <- mf
        coxph.call[[1L]] <- as.name("coxph")
        coxph.strata <- substitute(coxph.strata)
        coxph.call$data <- quote(coxph.data)
        coxph.data <- data
        if (!is.null(coxph.formula)) {
            coxph.formula2 <- coxph.call$formula
            rhs(coxph.formula2) <- rhs(formula) %call+% rhs(coxph.formula)
            coxph.call$formula <- coxph.formula2
        }
        if (!is.null(coxph.strata)) {
            coxph.formula2 <- coxph.call$formula
            rhs(coxph.formula2) <- rhs(formula) %call+% call("strata",coxph.strata)
            coxph.call$formula <- coxph.formula2
        }
        coxph.call$model <- TRUE
        ## coxph.obj <- eval(coxph.call, envir=parent.frame())
        coxph.obj <- eval(coxph.call, coxph.data)
        y <- model.extract(model.frame(coxph.obj),"response")
        data$logHhat <- pmax(-18,link$link(Shat(coxph.obj)))
    }
    if (interval) {
        ## survref regression
        survreg.call <- mf
        survreg.call[[1L]] <- as.name("survreg")
        survreg.obj <- eval(survreg.call, envir=parent.frame())
        weibullShape <- 1/survreg.obj$scale
        weibullScale <- predict(survreg.obj)
        y <- model.extract(model.frame(survreg.obj),"response")
        data$logHhat <- pmax(-18,link$link(pweibull(time,weibullShape,weibullScale,lower.tail=FALSE)))
    }
    ##
    ## initial values and object for lpmatrix predictions
    gam.call <- mf
    gam.call[[1L]] <- as.name("gam")
    gam.formula <- full.formula
    lhs(gam.formula) <- quote(logHhat) # new response
    gam.call$formula <- gam.formula
    gam.call$sp <- sp
    if (is.null(sp) && !is.null(sp.init) && (length(sp.init)>1 || sp.init!=1))
        gam.call$sp <- sp.init
    dataEvents <- data[event,]
    if (interval)
        dataEvents <- data
    gam.call$data <- quote(dataEvents) # events only
    gam.obj <- eval(gam.call)
    ## re-run gam if sp.init==1 (default)
    if (is.null(sp) && !is.null(sp.init) && length(sp.init)==1 && sp.init==1) {
        sp.init <- gam.call$sp <- rep(sp.init,length=length(gam.obj$sp))
        gam.obj <- eval(gam.call)
    }
    ##
    ## set up X, mf and wt
    mt <- terms(gam.obj)
    mf <- model.frame(gam.obj)
    wt <- if (is.null(substitute(weights))) rep(1,nrow(data)) else eval(substitute(weights),data,parent.frame())
    lpfunc <- function(x,...) {
      newdata <- data
      newdata[[timeVar]] <- x
      predict(gam.obj,newdata,type="lpmatrix")
    }
    ##
    bhazard <- substitute(bhazard)
    bhazard <- if (is.null(bhazard)) rep(0,nrow(data)) else eval(bhazard,data,parent.frame())
    excess <- !all(bhazard==0)
    ## initialise values specific to either delayed entry or interval-censored
    ind0 <- FALSE
    map0 <- 0L
    which0 <- 0
    wt0 <- 0
    ttype <- 0
    transX <- function(X, data) X
    transXD <- function(XD) XD
    ## browser()
    smooth <- gam.obj$smooth
    if (!interval) { # surv.type %in% c("right","counting")
        X <- predict(gam.obj,data,type="lpmatrix")
        if (link.type=="AH") {
            datat0 <- data
            datat0[[timeVar]] <- 0
            index0 <- which.dim(X - predict(gam.obj, datat0, type="lpmatrix"))
            smooth <- lapply(smooth, function(smoothi) {
                Sindex <- which((1:ncol(X) %in% index0)[smoothi$first.para:smoothi$last.para])
                para <- range(which((1:ncol(X) %in% smoothi$first.para:smoothi$last.para)[index0]))
                smoothi$S[[1]] <- smoothi$S[[1]][Sindex,Sindex]
                smoothi$first.para <- para[1]
                smoothi$last.para <- para[2]
                smoothi
                })
            transX <- function(X, data) {
                datat0 <- data
                datat0[[timeVar]] <- 0
                Xt0 <- predict(gam.obj, datat0, type="lpmatrix")
                (X - Xt0)[, index0, drop=FALSE]
            }
            transXD <- function(XD) XD[, index0, drop=FALSE]
            ## init <- init[index0]
        }
        X <- transX(X,data)
        XD <- grad1(lpfunc,data[[timeVar]])    
        XD <- transXD(matrix(XD,nrow=nrow(X)))
        X1 <- matrix(0,nrow(X),ncol(X))
        X0 <- matrix(0,1,ncol(X))
        if (delayed && all(time0==0)) delayed <- FALSE # CAREFUL HERE: delayed redefined
        if (delayed) {
            ind0 <- time0>0
            map0 <- vector("integer",nrow(X))
            map0[ind0] <- as.integer(1:sum(ind0))
            map0[!ind0] <- NaN
            ## which0 <- which(ind0)
            which0 <- 1:nrow(X)
            which0[!ind0] <- NaN
            data0 <- data[ind0,,drop=FALSE] # data for delayed entry times
            data0[[timeVar]] <- data0[[time0Var]]
            X0 <- transX(predict(gam.obj,data0,type="lpmatrix"), data0)
            wt0 <- wt[ind0]
            rm(data0)
        }
    } else { ## interval-censored
        ## ttime <- eventInstance[,1]
        ## ttime2 <- eventInstance[,2]
        ttype <- eventInstance[,3]
        X1 <- transX(predict(gam.obj,data,type="lpmatrix"), data)
        data0 <- data
        data0[[timeVar]] <- data0[[time0Var]]
        lpfunc <- function(x,...) {
            newdata <- data0
            newdata[[timeVar]] <- x
            predict(gam.obj,newdata,type="lpmatrix")
        }
        X <- transX(predict(gam.obj,data0,type="lpmatrix"), data0)
        XD <- grad1(lpfunc,data0[[timeVar]])    
        XD <- transXD(matrix(XD,nrow=nrow(X)))
        X0 <- matrix(0,nrow(X),ncol(X))
        rm(data0)
    }
    ## initial values
    if (is.null(init)) {
        init <- coef(gam.obj)
    }
    if (link.type=="AH") {
        init <- init[index0]
        }
    if (frailty) {
        init <- c(init,logtheta=logtheta)
    }
    ## smoothing parameters
    ## cases: 
    ##  (1) sp fixed
    ##  (2) sp.init
    ##  (3) use GAM
    if (no.sp <- is.null(sp)) {
        sp <- if(is.null(gam.obj$full.sp)) gam.obj$sp else gam.obj$full.sp
        if (!is.null(sp.init)) sp <- sp.init
    }
    if (!is.null(control) && "parscale" %in% names(control)) {
      if (length(control$parscale)==1)
        control$parscale <- rep(control$parscale,length(init))
      if (is.null(names(control$parscale)))
        names(control$parscale) <- names(init)
    } else {
        if(is.null(control)) 
            control <- list()
        control$parscale <- rep(1,length(init))
        names(control$parscale) <- names(init)
    }
    args <- list(init=init,X=X,XD=XD,bhazard=bhazard,wt=wt,event=ifelse(event,1,0),time=time,
                 delayed=delayed, interval=interval, X0=X0, wt0=wt0, X1=X1, parscale=control$parscale,
                 smooth=if(penalty == "logH") smooth else design,
                 sp=sp, reltol_search=reltol$search, reltol=reltol$final, reltol_outer=reltol$outer, trace=trace,
                 kappa=1.0,outer_optim=outer_optim,
                 alpha=alpha,criterion=switch(criterion,GCV=1,BIC=2),
                 oldcluster=cluster, cluster=if(!is.null(cluster)) as.vector(unclass(factor(cluster))) else NULL, frailty=frailty,
                 map0 = map0 - 1L, ind0 = ind0, which0=which0 - 1L, link = link.type,
                 penalty = penalty, ttype=ttype, RandDist=RandDist, optimiser=optimiser,
                 type=if (frailty && RandDist=="Gamma") "pstpm2_gamma_frailty" else if (frailty && RandDist=="LogN") "pstpm2_normal_frailty" else "pstpm2", recurrent = recurrent, maxkappa=maxkappa,
                 transX=transX, transXD=transXD, Z.formula = Z, thetaAO = theta.AO, excess=excess,
                 return_type="optim", data=data, robust_initial=robust_initial)
    if (frailty) {
        rule <- fastGHQuad::gaussHermiteData(nodes)
        args$gauss_x <- rule$x
        args$gauss_w <- rule$w
        args$adaptive <- adaptive
        args$Z <- model.matrix(Z, data)
        if (ncol(args$Z)>1) stop("Current implementation only allows for a single random effect")
        args$Z <- as.vector(args$Z)
    }
    ## penalty function
    pfun <- function(beta,sp) {
        sum(sapply(1:length(gam.obj$smooth),
                   function(i) {
                     smoother <- gam.obj$smooth[[i]]
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
        for (i in 1:length(gam.obj$smooth))
        {
          smoother <- gam.obj$smooth[[i]]
          ind <- smoother$first.para:smoother$last.para
          deriv[ind] <- sp[i] * smoother$S[[1]] %*% beta[ind]
        }
        return(deriv)
    }
    if (penalty == "h") {
        ## a current limitation is that the hazard penalty needs to extract the variable names from the smoother objects (e.g. log(time) will not work)
        stopifnot(sapply(gam.obj$smooth,function(obj) obj$term) %in% names(data) ||
                  !is.null(smoother.parameters))
        ## new penalty using the second derivative of the hazard
        design <- smootherDesign(gam.obj,data,smoother.parameters)
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
            if (frailty) beta <- beta[-length(beta)]
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
      if(!is.null(gam.obj$full.sp)) gam.obj$sp <- gam.obj$full.sp
      value <- NULL
      while(is.na(value <- negllsp(init,gam.obj$sp)) || !attr(value,"feasible")) {
          gam.call$sp <- gam.obj$sp * 5
          if (no.sp) sp <- gam.call$sp
          ## Unresolved: should we change sp.init if the initial values are not feasible?
          gam.obj <- eval(gam.call)
          if(!is.null(gam.obj$full.sp)) gam.obj$sp <- gam.obj$full.sp
          init <- coef(gam.obj)
          if (frailty)
              init <- c(init,logtheta=logtheta)
          if (all(gam.obj$sp > 1e5)) break
          ## stop("Initial values not valid and revised sp>1e5")
      }
      args$sp <- gam.obj$sp
    } else args$sp <- sp
#     ### Using exterior penalty method for nonlinear constraints: h(t)>=0 or increasing logH(t)
#     ### Some initial values should be outside the feasible region
#     while(all(XD%*%init>=0)){
#       init <- init+0.001
#     }
#     ### Check initial value
#     if(any(XD%*%init<=0)) {
#       cat("Some initial values are exactly outside the feasible region of this problem","\n") 
#     }
    ## MLE
    args$return_type <- if (!no.sp) { # fixed sp as specified
        args$return_type <- "optim_fixed"
    } else if (length(sp)>1) {
        "optim_multivariate"
    } else {
        "optim_first"
    }
    fit <- .Call("model_output", args, PACKAGE = "rstpm2")
    fit$coef <- as.vector(fit$coef)
    fit$sp <- as.vector(fit$sp)
    names(fit$coef) <- names(init)
    args$init <- init <- fit$coef
    args$sp <- sp <- fit$sp
    edf <- fit$edf
    edf_var<- as.vector(fit$edf_var)
    names(edf_var) <- sapply(gam.obj$smooth,"[[","label")
    names(fit$coef) <- rownames(fit$hessian) <- colnames(fit$hessian) <- names(init)
    args$kappa.final <- fit$kappa
    negll <- function(beta) negllsp(beta,sp)
    gradnegll <- function(beta) gradnegllsp(beta,sp)
    parnames(negll) <- parnames(gradnegll) <- names(init)
    mle2 <- if (use.gr) {
            mle2(negll,init,vecpar=TRUE, control=control, gr=gradnegll, eval.only=TRUE, ...)
        } else mle2(negll,init,vecpar=TRUE, control=control, eval.only=TRUE, ...)
    mle2@details$hessian <- fit$hessian
    ## mle2@vcov <- solve(optimHess(coef(mle2),negll,gradnegll))
    mle2@vcov <- solve(fit$hessian)
    mle2@details$convergence <- 0
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
               gam = gam.obj,
               timeVar = timeVar,
               time0Var = time0Var,
               timeExpr = timeExpr,
               time0Expr = time0Expr,
               like = like,
               fullformula = fullformula,
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
    # if (robust) # kludge
    #     out@vcov <- sandwich.stpm2(out, cluster=cluster)
    if (robust && !frailty) {
    	## Bread matrix
    	bread.mat <- solve(fit$hessian) 
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
    return(out)
}
## Could this inherit from summary.stpm2?
setClass("summary.pstpm2", representation(pstpm2="pstpm2",frailty="logical",theta="list",wald="matrix"), contains="summary.mle2")
setMethod("summary", "pstpm2",
          function(object) {
              newobj <- as(summary(as(object,"mle2")),"summary.pstpm2")
              newobj@pstpm2 <- object
              newobj@frailty <- object@frailty
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
                  statistic <- as.vector(coef1[i] %*% solve(vcov1[i,i]) %*% coef1[i])
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
          function(object,fun,newdata=NULL,link=c("I","log","cloglog","logit"),...)
  {
    link <- match.arg(link)
    invlinkf <- switch(link,I=I,log=exp,cloglog=cexpexp,logit=expit)
    linkf <- eval(parse(text=link))
    if (is.null(newdata) && !is.null(object@data))
      newdata <- object@data
    localf <- function(coef,...)
      {
        object@fullcoef = coef
        linkf(fun(object,...))
      }
    dm <- numDeltaMethod(object,localf,newdata=newdata,...)
    out <- invlinkf(data.frame(Estimate=dm$Estimate,
                               lower=dm$Estimate-1.96*dm$SE,
                               upper=dm$Estimate+1.96*dm$SE))
    ## cloglog switches the bounds
    if (link=="cloglog") 
      out <- data.frame(Estimate=out$Estimate,lower=out$upper,upper=out$lower)
    return(out)
  })
##
setMethod("predict", "pstpm2",
          function(object,newdata=NULL,
                   type=c("surv","cumhaz","hazard","density","hr","sdiff","hdiff","loghazard","link","meansurv","meansurvdiff","odds","or","margsurv","marghaz","marghr","meanhaz","af","fail","margfail","meanmargsurv"),
                   grid=FALSE,seqLength=300,
                   se.fit=FALSE,link=NULL,exposed=incrVar(var),var=NULL,keep.attributes=TRUE,use.gr=TRUE,...)
              predict.stpm2.base(object=object, newdata=newdata, type=type, grid=grid, seqLength=seqLength, se.fit=se.fit,
                                 link=link, exposed=exposed, var=var, keep.attributes=keep.attributes, use.gr=use.gr, ...))
setMethod("residuals", "pstpm2",
          function(object, type=c("li","gradli"))
              residuals.stpm2.base(object=object, type=type))

## setMethod("predict", "pstpm2",
##           function(object,newdata=NULL,
##                    type=c("surv","cumhaz","hazard","density","hr","sdiff","hdiff","loghazard","link","meansurv","meansurvdiff","odds","or","margsurv","marghaz","marghr","meanhaz"),
##                    grid=FALSE,seqLength=300,
##                    se.fit=FALSE,link=NULL,exposed=incrVar(var),var,...)
##   {
##     type <- match.arg(type)
##     ## exposed is a function that takes newdata and returns the revised newdata
##     ## var is a string for a variable that defines a unit change in exposure
##     local <-  function (object, newdata=NULL, type="surv", exposed)
##     {
##         args <- object@args
##         tt <- object@terms 
##         link <- object@link
##         if (is.null(newdata)) {
##           ##mm <- X <- model.matrix(object) # fails (missing timevar)
##           X <- object@x
##           XD <- object@xd
##           ##y <- model.response(object@model.frame)
##           y <- object@y
##           time <- as.vector(y[,ncol(y)-1])
##         }
##         else {
##           X <- args$transX(predict(object@gam, newdata, type="lpmatrix"), newdata)
##           ## lpfunc <- function(delta,fit,data,var) {
##           ##   data[[var]] <- data[[var]]+delta
##           ##   predict(fit,data,type="lpmatrix")
##           ## }
##           ## XD <- grad(lpfunc,0,object@gam,newdata,object@timeVar)
##           ## XD <- matrix(XD,nrow=nrow(X))
##           lpfunc <- function(x,...) {
##             newdata2 <- newdata
##             newdata2[[object@timeVar]] <- x
##             predict(object@gam,newdata2,type="lpmatrix")
##           }
##           XD <- args$transXD(grad1(lpfunc,newdata[[object@timeVar]]))
##           ## resp <- attr(Terms, "variables")[attr(Terms, "response")] 
##           ## similarly for the derivatives
##           ## if (object@delayed) {
##           ##   newdata0 <- newdata
##           ##   newdata0[[object@timeVar]] <- newdata[[object@time0Var]]
##           ##   X0 <- predict(object@gam, newdata0,type="lpmatrix")
##           ##   ## XD0 <- grad(lpfunc,0,object@lm,newdata,object@timeVar)
##           ##   ## XD0 <- matrix(XD0,nrow=nrow(X0))
##           ## }
##           if (type %in% c("hazard","hr","sdiff","hdiff","loghazard","meansurvdiff","or",
##                           "marghaz","marghr")) {
##             time <- eval(object@timeExpr,newdata)
##             ##
##           }
##           if (type %in% c("hr","sdiff","hdiff","meansurvdiff","or","marghr")) {
##             if (missing(exposed))
##               stop("exposed needs to be specified for type in ('hr','sdiff','hdiff','meansurvdiff','or','marghr')")
##             newdata2 <- exposed(newdata)
##             lpfunc <- function(x,...) {
##                 newdata3 <- newdata2
##                 newdata3[[object@timeVar]] <- x
##                 predict(object@gam,newdata3,type="lpmatrix")
##             }
##             X2 <- args$transX(predict(object@gam, newdata2, type="lpmatrix"), newdata2)
##             XD2 <- args$transXD(grad1(lpfunc,newdata2[[object@timeVar]]))
##             ## XD2 <- grad(lpfunc,0,object@gam,newdata2,object@timeVar)
##             ## XD2 <- matrix(XD2,nrow=nrow(X))
##           }
##         }
##         beta <- coef(object)
##         if (object@frailty) {
##             theta <- exp(beta[length(beta)])
##             beta <- beta[-length(beta)]
##             gauss_x <- object@args$gauss_x
##             gauss_w <- object@args$gauss_w
##             Z <- model.matrix(args$Z.formula, newdata)
##             if (ncol(Z)>1) stop("Current implementation only allows for a single random effect")
##             Z <- as.vector(Z)
##         }
##         eta <- as.vector(X %*% beta)
##         etaD <- as.vector(XD %*% beta)
##         S <- link$ilink(eta)
##         h <- link$h(eta,etaD)
##         if (any(h<0)) warning(sprintf("Predicted hazards less than zero (n=%i).",sum(h<0)))
##         H = link$H(eta)
##         Sigma = vcov(object)
##         if (type=="link") { # delayed entry?
##           return(eta)
##         }
##         if (type=="density")
##             return (S*h)
##         if (type=="cumhaz") { # delayed entry?
##             return(H)
##         }
##         if (type=="surv") { # delayed entry?
##           return(S)
##         }
##         if (type=="odds") { # delayed entry?
##           return((1-S)/S)
##         }
##         if (type=="sdiff")
##           return(link$ilink(as.vector(X2 %*% beta)) - S)
##         if (type=="or") {
##             S2 <- link$ilink(as.vector(X2 %*% beta)) 
##             return((1-S2)/S2/((1-S)/S))
##         }
##         if (type=="hazard") {
##           return(h)
##         }
##         if (type=="loghazard") {
##           return(log(h))
##         }
##         if (type=="hdiff") {
##             eta2 <- as.vector(X2 %*% beta)
##             etaD2 <- as.vector(XD2 %*% beta)
##             h2 <- link$h(eta2,etaD2)
##             return(h2 - h)
##         }
##         if (type=="hr") {
##             eta2 <- as.vector(X2 %*% beta)
##             etaD2 <- as.vector(XD2 %*% beta)
##             h2 <- link$h(eta2,etaD2)
##           return(h2/h)
##         }
##         if (type=="meansurv") {
##             return(mean(S))
##         }
##         if (type=="meanhaz") {
##             return(sum(h*S)/sum(S))
##         }
##         if (type=="meansurvdiff") {
##             eta2 <- as.vector(X2 %*% beta)
##             S2 <- link$ilink(eta2)
##             return(mean(S2-S))
##         }
##         if (type=="margsurv") {
##             stopifnot(object@frailty && object@args$RandDist %in% c("Gamma","LogN"))
##             if (object@args$RandDist=="Gamma")
##                 return((1+theta*H)^(-1/theta))
##             if (object@args$RandDist=="LogN") {
##                 return(sapply(1:length(gauss_x),
##                               function(i) link$ilink(eta+Z*sqrt(2)*sqrt(theta)*gauss_x[i])) %*%
##                        gauss_w / sqrt(pi))
##             }
##         }
##         if (type=="marghaz") {
##             stopifnot(object@frailty && object@args$RandDist %in% c("Gamma","LogN"))
##             if (object@args$RandDist=="Gamma") {
##                 margsurv <- (1+theta*H)^(-1/theta)
##                 return(h*margsurv^theta)
##             }
##             if (object@args$RandDist=="LogN") {
##                 return(sapply(1:length(gauss_x),
##                               function(i) link$h(eta+Z*sqrt(2)*sqrt(theta)*gauss_x[i],etaD)) %*%
##                        gauss_w / sqrt(pi))
##             }
##         }
##         if (type=="marghr") {
##             stopifnot(object@frailty && object@args$RandDist %in% c("Gamma","LogN"))
##             eta2 <- as.vector(X2 %*% beta)
##             etaD2 <- as.vector(XD2 %*% beta)
##             if (object@args$RandDist=="Gamma") {
##                 H2 <- link$H(eta2)
##                 h2 <- link$h(eta2,etaD2)
##                 margsurv <- (1+theta*H)^(-1/theta)
##                 marghaz <- h*margsurv^theta
##                 margsurv2 <- (1+theta*H2)^(-1/theta)
##                 marghaz2 <- h2*margsurv2^theta
##             }
##             if (object@args$RandDist=="LogN") {
##                 marghaz <- sapply(1:length(gauss_x),
##                                   function(i) as.vector(link$h(eta+Z*sqrt(2)*sqrt(theta)*gauss_x[i],etaD))) %*%
##                                   gauss_w / sqrt(pi)
##                 marghaz2 <- sapply(1:length(gauss_x),
##                                    function(i) as.vector(link$h(eta2+Z*sqrt(2)*sqrt(theta)*gauss_x[i],etaD2))) %*%
##                                        gauss_w / sqrt(pi)
##             }
##             return(marghaz2/marghaz)
##         }
##       }
##     ##debug(local)
##     if (is.null(newdata) && type %in% c("hr","sdiff","hdiff","meansurvdiff","or","marghr"))
##       stop("Prediction using type in ('hr','sdiff','hdiff','meansurvdiff','or','marghr') requires newdata to be specified.")
##     if (grid) {
##       Y <- object@y
##       event <- Y[,ncol(Y)]==1
##       time <- object@data[[object@timeVar]]
##       eventTimes <- time[event]
##       X <- seq(min(eventTimes),max(eventTimes),length=seqLength)[-1]
##       data.x <- data.frame(X)
##       names(data.x) <- object@timeVar
##       newdata <- merge(newdata,data.x)
##     }
##     pred <- if (!se.fit) {
##       local(object,newdata,type=type,exposed=exposed,
##             ...)
##     }
##     else {
##         if (is.null(link))
##         link <- switch(type,surv="cloglog",density="log",cumhaz="log",hazard="log",hr="log",sdiff="I",
##                        hdiff="I",loghazard="I",link="I",odds="log",or="log",margsurv="cloglog",marghaz="log",marghr="log",meanhaz="I",meansurv="I")
##       predictnl(object,local,link=link,newdata=newdata,type=type,
##                 exposed=exposed,...) 
##     }
##     attr(pred,"newdata") <- newdata
##     ##if (grid) cbind(newdata,as.data.frame(pred)) else pred
##     return(pred)
##   })

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
