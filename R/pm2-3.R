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
              centre=centre, log=log)
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
## first derivative of the nsx() function: deprecated (only used in rstpm2:::strpm2Old)
nsxDeriv<- 
  function (x, dorder=1, ...) {
    stopifnot(dorder %in% 1:2)
    basis <- nsx(x, ...)
    if (dorder==1) {
      h <- .Machine$double.eps^(1/3)*ifelse(abs(x)>1,abs(x),1)
      basisD <- apply(predict(basis,x+h)-predict(basis,x-h),2,
                      function(x) x/(2*h))
      if( attr(basis,"log"))
        basisD <- apply(basisD,2,function(y) y/exp(x))
    }
    if (dorder==2) {
      h <- .Machine$double.eps^(1/6)*ifelse(abs(x)>1,abs(x),1)
      basisD <- apply(predict(basis,x+h)-2*predict(basis,x)+predict(basis,x-h),2,
                      function(x) x/(h^2))
      if( attr(basis,"log")) {
        h <- .Machine$double.eps^(1/3)*ifelse(abs(x)>1,abs(x),1)
        basisD1 <- apply(predict(basis,x+h)-predict(basis,x-h),2,
                      function(x) x/(2*h))
        basisD <- apply(basisD-basisD1,2,function(y) y/exp(2*x))
      }
    }    
    attributes(basisD) <- attributes(basis)
    class(basisD) <- c("nsxDeriv", "basis")
    basisD
  }
## nsxDerivOld <- 
##   function (x, df = NULL, knots = NULL, intercept = FALSE,
##             Boundary.knots = range(x),
##             derivs = if (cure) c(2,1) else c(2,2),
##             centre = FALSE, cure = FALSE) {
##     mf <- match.call(expand.dots = FALSE)
##     mf[[1L]] <- as.name("nsx")
##     basis <- eval(mf, parent.frame())
##     h <- .Machine$double.eps^(1/3)*ifelse(abs(x)>1,abs(x),1)
##     basisD <- apply(predict(basis,x+h)-predict(basis,x-h),2,
##                       function(x) x/(2*h))
##     attributes(basisD) <- attributes(basis)
##     class(basisD) <- c("nsxDeriv", "basis")
##     basisD
##   }
## The following is deprecated (only used in rstpm2:::strpm2Old)
makepredictcall.nsxDeriv <- 
  function (var, call) 
{
  if (as.character(call)[1L] != "nsxDeriv") 
    return(call)
  at <- attributes(var)[c("knots", "Boundary.knots", "intercept",
                          "derivs", "centre", "log")]
  xxx <- call[1L:2]
  xxx[names(at)] <- at
  xxx
}
## The following is deprecated (only used in rstpm2:::strpm2Old)
predict.nsxDeriv <- 
  function (object, newx, ...) 
{
  if (missing(newx)) 
    return(object)
  a <- c(list(x = newx), attributes(object)[c("knots", "Boundary.knots", 
                "intercept", "derivs", "centre", "log")])
  do.call("nsxDeriv", a)
}
## nsDeriv<- 
## function (x, df = NULL, knots = NULL, intercept = FALSE, Boundary.knots = range(x)) 
## {
##     nx <- names(x)
##     x <- as.vector(x)
##     nax <- is.na(x)
##     if (nas <- any(nax)) 
##         x <- x[!nax]
##     if (!missing(Boundary.knots)) {
##         Boundary.knots <- sort(Boundary.knots)
##         outside <- (ol <- x < Boundary.knots[1L]) | (or <- x > 
##             Boundary.knots[2L])
##     }
##     else outside <- FALSE
##     if (!missing(df) && missing(knots)) {
##         nIknots <- df - 1 - intercept
##         if (nIknots < 0) {
##             nIknots <- 0
##             warning("'df' was too small; have used ", 1 + intercept)
##         }
##         knots <- if (nIknots > 0) {
##             knots <- seq.int(0, 1, length.out = nIknots + 2L)[-c(1L, 
##                 nIknots + 2L)]
##             stats::quantile(x[!outside], knots)
##         }
##     }
##     else nIknots <- length(knots)
##     Aknots <- sort(c(rep(Boundary.knots, 4L), knots))
##     if (any(outside)) {
##         basis <- array(0, c(length(x), nIknots + 4L))
##         basisD <- array(0, c(length(x), nIknots + 4L))
##         if (any(ol)) {
##             k.pivot <- Boundary.knots[1L]
##             xl <- cbind(1, x[ol] - k.pivot)
##             xlD <- cbind(0, rep(1,sum(ol)))
##             tt <- splines:::spline.des(Aknots, rep(k.pivot, 2L), 4, c(0, 
##                 1))$design
##              basis[ol, ] <- xl %*% tt
##             basisD[ol, ] <- xlD %*% tt
##         }
##         if (any(or)) {
##             k.pivot <- Boundary.knots[2L]
##             xr <- cbind(1, x[or] - k.pivot)
##             xrD <- cbind(0, rep(1,sum(or)))
##             tt <- splines:::spline.des(Aknots, rep(k.pivot, 2L), 4, c(0, 
##                 1))$design
##              basis[or, ] <- xr %*% tt
##             basisD[or, ] <- xrD %*% tt
##         }
##         if (any(inside <- !outside)) {
##              basis[inside, ] <- splines:::spline.des(Aknots, x[inside], 
##                  4)$design
##             basisD[inside, ] <- splines:::spline.des(Aknots, x[inside], 
##                 4, rep(1,sum(inside)))$design
##           }
##     }
##     else {
##       basis <- splines:::spline.des(Aknots, x, 4)$design
##       basisD <- splines:::spline.des(Aknots, x, 4, rep(1,length(x)))$design
##     }
##     const <- splines:::spline.des(Aknots, Boundary.knots, 4, c(2, 2))$design
##     if (!intercept) {
##         const <- const[, -1, drop = FALSE]
##         basis <- basis[, -1, drop = FALSE]
##         basisD <- basisD[, -1, drop = FALSE]
##     }
##     qr.const <- qr(t(const))
##     basis <- as.matrix((t(qr.qty(qr.const, t(basis))))[, -(1L:2L), 
##         drop = FALSE])
##     basisD <- as.matrix((t(qr.qty(qr.const, t(basisD))))[, -(1L:2L), 
##         drop = FALSE])
##     n.col <- ncol(basis)
##     if (nas) {
##         nmat <- matrix(NA, length(nax), n.col)
##         nmat[!nax, ] <- basisD
##         basisD <- nmat
##     }
##     dimnames(basisD) <- list(nx, 1L:n.col)
##     a <- list(degree = 3, knots = if (is.null(knots)) numeric(0) else knots, 
##         Boundary.knots = Boundary.knots, intercept = intercept)
##     attributes(basisD) <- c(attributes(basisD), a)
##     class(basisD) <- c("nsDeriv", "basis")
##     basisD
## }
Shat <- function(obj)
  {
    ## predicted survival for individuals (adjusted for covariates)
    newobj = survfit(obj,se.fit=FALSE)
    surv = newobj$surv
    rr = try(predict(obj,type="risk"),silent=TRUE)
    ## case: only an intercept in the main formula with strata (it would be better to recognise this using attributes for newobj)
    if (inherits(rr,"try-error")) {
      if (try %in% c("Error in colSums(x[j, ] * weights[j]) : \n  'x' must be an array of at least two dimensions\n",
                  "Error in rowsum.default(x * weights, indx) : incorrect length for 'group'\n")) rr <- 1 else stop(rr)
    }
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
numDeltaMethod <- function(object,fun,...) {
  coef <- coef(object)
  est <- fun(coef,...)
  Sigma <- vcov(object)
  gd <- grad(fun,coef,...)
  se.est <- as.vector(sqrt(diag(t(gd) %*% Sigma %*% gd)))
  data.frame(Estimate = est, SE = se.est)
}
predictnl <- function (object, ...) 
  UseMethod("predictnl")
predictnl.default <- function(object,fun,newdata=NULL,...)
  {
    ## link=c(I,log,sqrt),invlink=NULL
    ## link <- match.arg(link)
    ## if (is.null(invlink))
    ##       invlink <- switch(deparse(substitute(link)),I=I,log=exp,sqrt=function(x) x^2)
    if (is.null(newdata) && !is.null(object$data))
      newdata <- object$data
    localf <- function(coef,...)
      {
        object$coefficients = coef
        fun(object,...)
      }
    numDeltaMethod(object,localf,newdata=newdata,...)
  }
setMethod("predictnl", "mle2", function(object,fun,newdata=NULL,...)
  {
    if (is.null(newdata) && !is.null(object@data))
      newdata <- object@data
    localf <- function(coef,...)
      {
        object@fullcoef = coef # changed from predictnl.default()
        fun(object,...)
      }
    numDeltaMethod(object,localf,newdata=newdata,...)
  })
## setMethod("predictnl", "mle", function(object,fun,...)
##   {
##     localf <- function(coef,...)
##       {
##         object@fullcoef = coef # changed from predictnl.default()
##         fun(object,...)
##       }
##     numDeltaMethod(object,localf,...)
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
  rr <- t(rstpm2:::grad(obj@logli,coef(obj)))
  colnames(rr) <- names(coef(obj))
  if (weighted)
    rr <- rr * obj@weights
  rr
}
meat.stpm2 <- 
function (x, adjust = FALSE, ...) 
{
    psi <- estfun.stpm2(x, ...)
    k <- NCOL(psi)
    n <- NROW(psi)
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
  if (n>1)
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
##setClassUnion("numericOrNULL",c("numeric","NULL"))
setOldClass("Surv")
setOldClass("lm")

## old stpm2 implementation (included for comparison purposes - deprecated)
setClass("stpm2Old", representation(xlevels="list",
                                 contrasts="listOrNULL",
                                 terms="terms",
                                 logli="function",
                                 ## weights="numericOrNULL",
                                 model.frame="list",
                                 x="matrix",
                                 xd="matrix",
                                 termsd="terms",
                                 Call="call",
                                 y="Surv"
                                 ),
         contains="mle2")
stpm2Old <- function(formula, data,
                  df = 3, cure = FALSE, logH.args = NULL, logH.formula = NULL,
                  tvc = NULL, tvc.formula = NULL,
                  control = list(parscale = 0.1, maxit = 300), init = FALSE,
                  coxph.strata = NULL, weights = NULL, robust = FALSE, baseoff = FALSE,
                  bhazard = NULL, contrasts = NULL, subset = NULL, ...)
  {
    ## ensure that data is a data frame
    data <- get_all_vars(formula, data)
    ## restrict to non-missing data (assumes na.action=na.omit)
    .include <- Reduce(`&`,
                       lapply(model.frame(formula, data, na.action=na.pass),
                              Negate(is.na)),
                       TRUE)
    data <- data[.include, , drop=FALSE]
    Call <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "contrasts", "weights"),
               names(mf), 0L)
    mf <- mf[c(1L, m)]
    ## mf$drop.unused.levels <- TRUE # include?
    mf[[1L]] <- as.name("model.frame")
    eventExpression <- lhs(formula)[[length(lhs(formula))]]
    delayed <- length(lhs(formula))==4
    timeExpr <- lhs(formula)[[if (delayed) 3 else 2]] # expression
    timeVar <- all.vars(timeExpr)
    stopifnot(length(timeVar)==1)
    ## set up the formulae
    if (is.null(logH.formula) && is.null(logH.args)) {
      logH.args$df <- df
      if (cure) logH.args$cure <- cure
    }
    if (!is.null(logH.args) && is.null(logH.args$log))
      logH.args$log <- TRUE
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
                              vector2call(list(log=TRUE,cure=cure,df=tvc[[name]]))))))
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
    ## differentiation would be easier if we did not have functions for log(time)
    logHD.formula <- replaceFormula(logH.formula,quote(nsx),quote(nsxDeriv))
    ## set up primary terms objects (mt and mtd)
    mf$formula = full.formula
    datae <- data[eval(eventExpression,data)==1, , drop=FALSE]
    mf$data <- quote(datae) # restricted to event times
    mfX <- mfd <- mf # copy
    mf <- eval(mf)
    mt <- attr(mf, "terms") # primary!
    xlev <- .getXlevels(mt, mf)
    mfd[[2]] <- logHD.formula
    mfd <- eval(mfd)
    mtd <- attr(mfd, "terms") # primary!
    ## design matrices
    mfX$formula <- quote(mt)
    mfX$data <- quote(data)
    mfX2 <- mfXD <- mfX # copies
    mfX <- eval(mfX)
    if (!is.null(cl <- attr(mt, "dataClasses"))) 
      .checkMFClasses(cl, mfX)
    X <- model.matrix(mt, mfX, contrasts)
    wt <- model.weights(mfX)
    if (is.null(wt)) wt <- rep(1,nrow(X))
    ## mfXD <- model.frame(mtd, data, xlev=xlev, weights = weights)
    mfXD$formula <- quote(mtd)
    mfXD <- eval(mfXD)
    if (!is.null(cl <- attr(mtd, "dataClasses"))) 
      .checkMFClasses(cl, mfXD)
    XD <- model.matrix(mtd, mfXD, contrasts)[,-1,drop=FALSE]
    ##
    y <- model.extract(mfX,"response")
    if (!inherits(y, "Surv")) 
      stop("Response must be a survival object")
    type <- attr(y, "type")
    if (type != "right" && type != "counting") 
        stop(paste("stpm2Old model doesn't support \"", type, "\" survival data", 
            sep = ""))
    event <- y[,ncol(y)]==1
    time <- y[,ncol(y)-1]
    ## initial values
    coxph.strata <- substitute(coxph.strata)
    coxph.formula <- formula
    if (!is.null(coxph.strata)) 
      rhs(coxph.formula) <- rhs(formula) %call+% call("strata",coxph.strata)
    coxph.obj <- quote(coxph(coxph.formula,data=data,model=TRUE))
    coxph.obj$weights <- substitute(weights)
    coxph.obj <- eval(coxph.obj)
    ## coxph.obj <- eval.parent(substitute(coxph(formula,data),
    ##                                     list(formula=formula,data=data)))
    ##coxph.data <- get_all_vars(mfX)
    coxph.data <- data
    coxph.data$logHhat <- pmax(-18,log(-log(Shat(coxph.obj))))
    coxph.data[[timeVar]] <- data[[timeVar]]
    ## coxph.data[,1] <- time
    ## names(coxph.data)[1] <- as.character(timeExpr)
    ##
    lm.formula <- full.formula
    lhs(lm.formula) <- quote(logHhat) # new response
    if (!init) {
      lm.obj <- lm(lm.formula,coxph.data[event,],contrasts=contrasts)
      ## lm.obj <- lm(lm.formula,data[event,],contrasts=contrasts)
      lm.obj$weights <- substitute(weights)
      lm.obj <- eval(lm.obj)
      init <- coef(lm.obj)
    }
    ## indexXD <- grep(timeVar,names(init)) ## FRAUGHT: pattern-matching on names
    data2 <- data
    data2[[timeVar]] <- data[[timeVar]]+1e-3
    temp <- predict(full.formula,data,data)
    temp2 <- predict(full.formula,data,data2)
    indexXD <- apply(apply(temp2-temp,2,range),2,diff)>1e-8
    rm(data2,temp,temp2)
    bhazard <- if (is.null(bhazard)) 0 else bhazard[event] # crude
    if (delayed && any(y[,1]>0)) {
      data2 <- data[y[,1]>0,,drop=FALSE] # data for delayed entry times
      mt2 <- delete.response(mt)
      ## copy over times - so we can use the same term object
      timeExpr0 <- lhs(full.formula)[[2]] # expression
      timeVar0 <- all.vars(timeExpr0)
      stopifnot(length(timeVar0)==1)
      data2[[timeVar]] <- data2[[timeVar0]] 
      ##y.names <- sapply(lhs(full.formula),deparse)
      ##data2[[y.names[3]]] <- data2[[y.names[2]]] 
      ## data2[[y.names[2]]] <- 0
      mfX2$formula <- quote(mt2)
      mfX2$data <- quote(data2)
      mfX2 <- eval(mfX2)
      ##mfX2 <- model.frame(mt2, data2, xlev=xlev, weights = weights)
      if (!is.null(cl <- attr(mt2, "dataClasses"))) 
        .checkMFClasses(cl, mfX2)
      X2 <- model.matrix(mt2, mfX2, contrasts)
      wt2 <- model.weights(mfX2)
      if (is.null(wt2)) wt2 <- rep(1,nrow(X2))
      ## delayed.formula <- full.formula
      ## lhs(delayed.formula) <- NULL
      ## X2 = predict(delayed.formula, data, data2) ## delayed entry
      negll <- function(beta) {
        eta <- X %*% beta
        eta2 <- X2 %*% beta
        h <- (XD[event,] %*% beta[indexXD])*exp(eta[event]) + bhazard
        ## h <- (XD[event,] %*% beta[indexXD])*exp(eta[event])/time[event] + bhazard
        h[h<0] <- 1e-100
        ll <- sum(wt[event]*log(h)) +  sum(wt2*exp(eta2)) -
          sum(wt*exp(eta))
        return(-ll)
      }
      logli <- function(beta) {
        eta <- X %*% beta
        eta2 <- X2 %*% beta
        h <- (XD[event,] %*% beta[indexXD])*exp(eta[event]) + bhazard
        h[h<0] <- 1e-100
        out <- exp(eta2) - exp(eta)
        out[event] <- out[event]+log(h)
        out <- out*wt
        return(out)
      }
    }
    else { # right censoring only
      negll <- function(beta) {
        eta <- X %*% beta
        h <- (XD[event,] %*% beta[indexXD])*exp(eta[event]) + bhazard
        ## h <- (XD[event,] %*% beta[indexXD])*exp(eta[event])/time[event] + bhazard
        h[h<0] <- 1e-100
        ll <- sum(wt[event]*log(h)) - sum(wt*exp(eta))
        return(-ll)
      }
      logli <- function(beta) {
        eta <- X %*% beta
        h <- (XD[event,] %*% beta[indexXD])*exp(eta[event]) + bhazard
        h[h<0] <- 1e-100
        out <- -exp(eta)
        out[event] <- out[event]+log(h)
        out <- out*wt
        return(out)
      }
    }
    ## MLE
    if (!is.null(control) && "parscale" %in% names(control)) {
      if (length(control$parscale)==1)
        control$parscale <- rep(control$parscale,length(init))
      if (is.null(names(control$parscale)))
        names(control$parscale) <- names(init)
    }
    parnames(negll) <- names(init)
    mle2 <- mle2(negll,init,vecpar=TRUE, control=control, ...)
    out <- new("stpm2Old",
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
               contrasts = attr(X, "contrasts"),
               logli = logli,
               ##weights = weights,
               Call = Call,
               terms = mt,
               model.frame = mf,
               x = X,
               xd = XD,
               termsd = mtd,
               y = y)
    if (robust) # kludge
      out@vcov <- sandwich.stpm2(out)
    return(out)
  }
##
## predict survival
setMethod("predictnl", "stpm2Old",
          function(object,fun,newdata=NULL,link=c("I","log","cloglog"),...)
  {
    link <- match.arg(link)
    invlinkf <- switch(link,I=I,log=exp,cloglog=cexpexp)
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
setMethod("predict", "stpm2Old",
          function(object,newdata=NULL,
                   type=c("surv","cumhaz","hazard","hr","sdiff","hdiff","loghazard","link"),
                   grid=FALSE,seqLength=300,
                   se.fit=FALSE,link=NULL,exposed=incrVar(var),var,...)
  {
    ## exposed is a function that takes newdata and returns the revised newdata
    ## var is a string for a variable that defines a unit change in exposure
    local <-  function (object, newdata=NULL, type="surv", exposed)
      {
        tt <- object@terms 
        if (is.null(newdata)) {
          ##mm <- X <- model.matrix(object) # fails (missing timevar)
          X <- object@x
          XD <- object@xd
          ##y <- model.response(object@model.frame)
          y <- object@y
          time <- as.vector(y[,ncol(y)-1])
        }
        else {
          Terms <- delete.response(tt)
          m <- model.frame(Terms, newdata, xlev = object@xlevels)
          if (!is.null(cl <- attr(Terms, "dataClasses"))) 
            .checkMFClasses(cl, m)
          X <- model.matrix(Terms, m, contrasts.arg = object@contrasts)
          resp <- attr(Terms, "variables")[attr(Terms, "response")]
          ## similarly for the derivatives
          if (type %in% c("hazard","hr","sdiff","hdiff","loghazard")) {
            ttd <- object@termsd
            TermsD <- delete.response(ttd)
            md <- model.frame(TermsD, newdata, xlev = object@xlevels)
            if (!is.null(cld <- attr(TermsD, "dataClasses"))) 
              .checkMFClasses(cld, md)
            XD <- model.matrix(TermsD, md, contrasts.arg = object@contrasts)[,-1,drop=FALSE]
            ## how to elegantly extract the time variable?
            timevar <- if (length(tt[[2]])==3) tt[[2]][[2]] else tt[[2]][[3]]
            time <- model.matrix(as.formula(call("~",timevar)),newdata)[,-1,drop=TRUE]
            ##
          }
          if (type %in% c("hr","sdiff","hdiff")) {
            if (missing(exposed))
              stop("exposed needs to be specified for type in ('hr','sdiff','hdiff')")
            newdata2 <- exposed(newdata)
            m2 <- model.frame(Terms, newdata2, xlev = object@xlevels)
            if (!is.null(cl <- attr(Terms, "dataClasses"))) 
              .checkMFClasses(cl, m2)
            X2 <- model.matrix(Terms, m2, contrasts.arg = object@contrasts)
            md2 <- model.frame(TermsD, newdata2, xlev = object@xlevels)
            if (!is.null(cld <- attr(TermsD, "dataClasses"))) 
              .checkMFClasses(cld, md2)
            XD2 <- model.matrix(TermsD, md2, contrasts.arg = object@contrasts)[,-1,drop=FALSE]
          }
        }
        beta <- coef(object)
        cumHaz = exp(X %*% beta)
        Sigma = vcov(object)
        if (type=="link") { # delayed entry?
          return(X %*% beta)
        }
        if (type=="cumhaz") { # delayed entry?
          return(cumHaz)
        }
        if (type=="surv") { # delayed entry?
          return(exp(-cumHaz))
        }
        if (type=="sdiff")
          return(exp(-exp(X2 %*% beta)) - exp(-cumHaz))
        if (type=="hazard") {
          betaXD <- beta[(ncol(X)-ncol(XD)+1):ncol(X)]
          ##return((XD %*% betaXD)/time*cumHaz)
          return((XD %*% betaXD)*cumHaz)
        }
        if (type=="loghazard") {
          betaXD <- beta[(ncol(X)-ncol(XD)+1):ncol(X)]
          ##return((XD %*% betaXD)/time*cumHaz)
          return(log(XD %*% betaXD)+log(cumHaz))
        }
        if (type=="hdiff") {
          betaXD <- beta[(ncol(X)-ncol(XD)+1):ncol(X)]
          ##return((XD2 %*% betaXD)/time*exp(X2 %*% beta) - (XD %*% betaXD)/time*cumHaz)
          return((XD2 %*% betaXD)*exp(X2 %*% beta) - (XD %*% betaXD)/time*cumHaz)
        }
        if (type=="hr") {
          cumHazRatio = exp((X2 - X) %*% beta)
          betaXD <- beta[(ncol(X)-ncol(XD)+1):ncol(X)]
          return((XD2 %*% betaXD)/(XD %*% betaXD)*cumHazRatio)
        }
      }
    ##debug(local)
    type <- match.arg(type)
    if (is.null(newdata) && type %in% c("hr","sdiff","hdiff"))
      stop("Prediction using type in ('hr','sdiff','hdiff') requires newdata to be specified.")
    if (grid) {
      Terms <- object@terms
      timevar <- lhs(Terms)[[length(lhs(Terms))-1]]
      Y <- object@y
      event <- Y[,ncol(Y)]==1
      time <- Y[,ncol(Y)-1]
      eventTimes <- time[event]
      X <- seq(min(eventTimes),max(eventTimes),length=seqLength)[-1]
      data.x <- data.frame(X)
      names(data.x) <- deparse(timevar)
      newdata <- merge(newdata,data.x)
    }
    pred <- if (!se.fit) {
      local(object,newdata,type=type,exposed=exposed,
            ...)
    }
    else {
      if (is.null(link))
        link <- switch(type,surv="cloglog",cumhaz="log",hazard="log",hr="log",sdiff="I",
                       hdiff="I",loghazard="I",link="I")
      predictnl(object,local,link=link,newdata=newdata,type=type,
                exposed=exposed,...) 
    }
    attr(pred,"newdata") <- newdata
    ##if (grid) cbind(newdata,as.data.frame(pred)) else pred
    return(pred)
  })
setMethod("plot", signature(x="stpm2Old", y="missing"),
          function(x,y,newdata,type="surv",
                      xlab="Time",ylab=NULL,line.col=1,ci.col="grey",lty=par("lty"),
                      add=FALSE,ci=TRUE,rug=TRUE,
                      var=NULL,...) {
  y <- predict(x,newdata,type=type,var=var,grid=TRUE,se.fit=TRUE)
  if (is.null(ylab))
    ylab <- switch(type,hr="Hazard ratio",hazard="Hazard",surv="Survival",
                   sdiff="Survival difference",hdiff="Hazard difference",cumhaz="Cumulative hazard")
  xx <- attr(y,"newdata")
  xx <- xx[,ncol(xx)]
  if (!add) matplot(xx, y, type="n", xlab=xlab, ylab=ylab, ...)
  if (ci) polygon(c(xx,rev(xx)), c(y[,2],rev(y[,3])), col=ci.col, border=ci.col)
  lines(xx,y[,1],col=line.col,lty=lty)
  if (rug) {
      Y <- x@y
      eventTimes <- Y[Y[,ncol(Y)]==1,ncol(Y)-1]
      rug(eventTimes,col=line.col)
    }
  return(invisible(y))
})

## re-factored stpm2
setClass("stpm2", representation(xlevels="list",
                                 contrasts="listOrNULL",
                                 terms="terms",
                                 logli="function",
                                 ## weights="numericOrNULL",
                                    lm="lm",
                                    timeVar="character",
                                    timeExpr="nameOrcall",
                                 model.frame="list",
                                    call.formula="formula",
                                 x="matrix",
                                 xd="matrix",
                                 termsd="terms",
                                 Call="call",
                                 y="Surv"
                                 ),
         contains="mle2")
## TODO: weights
stpm2 <- function(formula, data,
                     df = 3, cure = FALSE, logH.args = NULL, logH.formula = NULL,
                     tvc = NULL, tvc.formula = NULL,
                     control = list(parscale = 0.1, maxit = 300), init = FALSE,
                     coxph.strata = NULL, weights = NULL, robust = FALSE, baseoff = FALSE,
                     bhazard = NULL, timeVar = NULL, use.gr = TRUE,
                     contrasts = NULL, subset = NULL, ...)
  {
    ## set up the data
    ## ensure that data is a data frame
    data <- get_all_vars(formula, data)
    ## restrict to non-missing data (assumes na.action=na.omit)
    .include <- Reduce(`&`,
                       lapply(model.frame(formula, data, na.action=na.pass),
                              Negate(is.na)),
                       TRUE)
    data <- data[.include, , drop=FALSE]
    ##
    ## parse the function call
    Call <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "contrasts", "weights"),
               names(mf), 0L)
    mf <- mf[c(1L, m)]
    ##
    ## parse the event expression
    stopifnot(length(lhs(formula))>=2)
    eventExpr <- lhs(formula)[[length(lhs(formula))]]
    delayed <- length(lhs(formula))==4
    timeExpr <- lhs(formula)[[if (delayed) 3 else 2]] # expression
    if (is.null(timeVar))
      timeVar <- all.vars(timeExpr)
    time <- eval(timeExpr, data)
    if (delayed) {
      time0Expr <- lhs(formula)[[2]]
      time0 <- eval(time0Expr, data)
    }
    event <- eval(eventExpr,data)
    event <- event > min(event)
    ##
    ## set up the formulae
    if (is.null(logH.formula) && is.null(logH.args)) {
      logH.args$df <- df
      if (cure) logH.args$cure <- cure
    }
    if (!is.null(logH.args) && is.null(logH.args$log))
      logH.args$log <- TRUE
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
                              vector2call(list(log=TRUE,cure=cure,df=tvc[[name]]))))))
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
    ##
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
    data$logHhat <- pmax(-18,log(-log(Shat(coxph.obj))))
    ##
    ## initial values and object for lpmatrix predictions
    lm.call <- mf
    lm.call[[1L]] <- as.name("lm")
    lm.formula <- full.formula
    lhs(lm.formula) <- quote(logHhat) # new response
    lm.call$formula <- lm.formula
    dataEvents <- data[event,]
    lm.call$data <- quote(dataEvents) # events only
    lm.obj <- eval(lm.call)
    if (!init) {
      init <- coef(lm.obj)
    }
    ##
    ## set up X, mf and wt
    X <- lpmatrix.lm(lm.obj,data)
    mt <- terms(lm.obj)
    mf <- model.frame(lm.obj)
    wt <- model.weights(lm.obj$model)
    if (is.null(wt)) wt <- rep(1,nrow(X))
    ##
    ## XD matrix
    lpfunc <- function(delta,fit,data,var) {
      data[[var]] <- data[[var]]+delta
      lpmatrix.lm(fit,data)
    }
    XD <- grad(lpfunc,0,lm.obj,data,timeVar)
    XD <- matrix(XD,nrow=nrow(X))
    ##
    bhazard <- substitute(bhazard)
    bhazard <- if (is.null(bhazard)) 0 else eval(bhazard,data,parent.frame())
    if (delayed && any(time0>0)) {
      ind <- time0>0
      data2 <- data[ind,,drop=FALSE] # data for delayed entry times
      X2 <- lpmatrix.lm(lm.obj, data2)
      wt2 <- wt[ind]
      gradnegll <- function(beta) {
        eta <- as.vector(X %*% beta)
        eta2 <- as.vector(X2 %*% beta)
        etaD <- as.vector(XD %*% beta)
        h <- etaD*exp(eta) + bhazard
        h[h<0] <- 1e-100
        g <- colSums(exp(eta)*wt*(-X + ifelse(event,1/h,0)*(XD + X*etaD)))-
          colSums(exp(eta2)*wt2*X2)
        return(-g)
      }
      negll <- function(beta) {
        eta <- X %*% beta
        eta2 <- X2 %*% beta
        h <- (XD %*% beta)*exp(eta) + bhazard
        h[h<0] <- 1e-100
        ll <- sum(wt[event]*log(h[event])) +  sum(wt2*exp(eta2)) -
          sum(wt*exp(eta))
        return(-ll)
      }
      logli <- function(beta) {
        eta <- X %*% beta
        eta2 <- X2 %*% beta
        h <- (XD %*% beta)*exp(eta) + bhazard
        h[h<0] <- 1e-100
        out <- exp(eta2) - exp(eta)
        out[event] <- out[event]+log(h[event])
        out <- out*wt
        return(out)
      }
    }
    else { # right censoring only
      gradnegll <- function(beta) {
        eta <- as.vector(X %*% beta)
        etaD <- as.vector(XD %*% beta)
        h <- etaD*exp(eta) + bhazard
        h[h<0] <- 1e-100
        g <- colSums(exp(eta)*wt*(-X + ifelse(event,1/h,0)*(XD + X*etaD)))
        return(-g)
      }
      negll <- function(beta) {
        eta <- X %*% beta
        h <- (XD %*% beta)*exp(eta) + bhazard
        h[h<0] <- 1e-100
        ll <- sum(wt[event]*log(h[event])) - sum(wt*exp(eta))
        return(-ll)
      }
      logli <- function(beta) {
        eta <- X %*% beta
        h <- (XD %*% beta)*exp(eta) + bhazard
        h[h<0] <- 1e-100
        out <- -exp(eta)
        out[event] <- out[event]+log(h[event])
        out <- out*wt
        return(out)
      }
    }
    ## MLE
    if (!is.null(control) && "parscale" %in% names(control)) {
      if (length(control$parscale)==1)
        control$parscale <- rep(control$parscale,length(init))
      if (is.null(names(control$parscale)))
        names(control$parscale) <- names(init)
    }
    parnames(negll) <- parnames(gradnegll) <- names(init)
    mle2 <- if (use.gr)
      mle2(negll,init,vecpar=TRUE, control=control, gr=gradnegll, ...)
    else mle2(negll,init,vecpar=TRUE, control=control, ...)
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
               contrasts = NULL, # wrong!
               logli = logli,
               ##weights = weights,
               Call = Call,
               terms = mt,
               model.frame = mf,
               lm = lm.obj,
               timeVar = timeVar,
               timeExpr = timeExpr,
               call.formula = formula,
               x = X,
               xd = XD,
               termsd = mt, # wrong!
               y = y)
    if (robust) # kludge
      out@vcov <- sandwich.stpm2(out)
    return(out)
  }
setMethod("predictnl", "stpm2",
          function(object,fun,newdata=NULL,link=c("I","log","cloglog"),...)
  {
    link <- match.arg(link)
    invlinkf <- switch(link,I=I,log=exp,cloglog=cexpexp)
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
setMethod("predict", "stpm2",
          function(object,newdata=NULL,
                   type=c("surv","cumhaz","hazard","hr","sdiff","hdiff","loghazard","link"),
                   grid=FALSE,seqLength=300,
                   se.fit=FALSE,link=NULL,exposed=incrVar(var),var,...)
  {
    ## exposed is a function that takes newdata and returns the revised newdata
    ## var is a string for a variable that defines a unit change in exposure
    local <-  function (object, newdata=NULL, type="surv", exposed)
      {
        tt <- object@terms 
        if (is.null(newdata)) {
          ##mm <- X <- model.matrix(object) # fails (missing timevar)
          X <- object@x
          XD <- object@xd
          ##y <- model.response(object@model.frame)
          y <- object@y
          time <- as.vector(y[,ncol(y)-1])
        }
        else {
          lpfunc <- function(delta,fit,data,var) {
            data[[var]] <- data[[var]]+delta
            lpmatrix.lm(fit,data)
          }
          X <- lpmatrix.lm(object@lm, newdata)
          XD <- grad(lpfunc,0,object@lm,newdata,object@timeVar)
          XD <- matrix(XD,nrow=nrow(X))
          ## resp <- attr(Terms, "variables")[attr(Terms, "response")] 
          ## similarly for the derivatives
          if (type %in% c("hazard","hr","sdiff","hdiff","loghazard")) {
            ## how to elegantly extract the time variable?
            timeExpr <- 
              lhs(object@call.formula)[[length(lhs(object@call.formula))-1]]
            time <- eval(timeExpr,newdata)
            ##
          }
          if (type %in% c("hr","sdiff","hdiff")) {
            if (missing(exposed))
              stop("exposed needs to be specified for type in ('hr','sdiff','hdiff')")
            newdata2 <- exposed(newdata)
            X2 <- lpmatrix.lm(object@lm, newdata2)
            XD2 <- grad(lpfunc,0,object@lm,newdata2,object@timeVar)
            XD2 <- matrix(XD,nrow=nrow(X))
          }
        }
        beta <- coef(object)
        cumHaz = exp(X %*% beta)
        Sigma = vcov(object)
        if (type=="link") { # delayed entry?
          return(X %*% beta)
        }
        if (type=="cumhaz") { # delayed entry?
          return(cumHaz)
        }
        if (type=="surv") { # delayed entry?
          return(exp(-cumHaz))
        }
        if (type=="sdiff")
          return(exp(-exp(X2 %*% beta)) - exp(-cumHaz))
        if (type=="hazard") {
          return((XD %*% beta)*cumHaz)
        }
        if (type=="loghazard") {
          return(log(XD %*% beta)+log(cumHaz))
        }
        if (type=="hdiff") {
          return((XD2 %*% beta)*exp(X2 %*% beta) - (XD %*% betaXD)/time*cumHaz)
        }
        if (type=="hr") {
          cumHazRatio = exp((X2 - X) %*% beta)
          return((XD2 %*% beta)/(XD %*% beta)*cumHazRatio)
        }
      }
    ##debug(local)
    type <- match.arg(type)
    if (is.null(newdata) && type %in% c("hr","sdiff","hdiff"))
      stop("Prediction using type in ('hr','sdiff','hdiff') requires newdata to be specified.")
    if (grid) {
      Y <- object@y
      event <- Y[,ncol(Y)]==1
      time <- object@data[[object@timeVar]]
      eventTimes <- time[event]
      X <- seq(min(eventTimes),max(eventTimes),length=seqLength)[-1]
      data.x <- data.frame(X)
      names(data.x) <- object@timeVar
      newdata <- merge(newdata,data.x)
    }
    pred <- if (!se.fit) {
      local(object,newdata,type=type,exposed=exposed,
            ...)
    }
    else {
      if (is.null(link))
        link <- switch(type,surv="cloglog",cumhaz="log",hazard="log",hr="log",sdiff="I",
                       hdiff="I",loghazard="I",link="I")
      predictnl(object,local,link=link,newdata=newdata,type=type,
                exposed=exposed,...) 
    }
    attr(pred,"newdata") <- newdata
    ##if (grid) cbind(newdata,as.data.frame(pred)) else pred
    return(pred)
  })
##`%c%` <- function(f,g) function(...) g(f(...)) # function composition
setMethod("plot", signature(x="stpm2", y="missing"),
          function(x,y,newdata,type="surv",
                      xlab=NULL,ylab=NULL,line.col=1,ci.col="grey",lty=par("lty"),
                      add=FALSE,ci=TRUE,rug=TRUE,
                      var=NULL,...) {
  y <- predict(x,newdata,type=type,var=var,grid=TRUE,se.fit=TRUE)
  if (is.null(xlab)) xlab <- deparse(x@timeExpr)
  if (is.null(ylab))
    ylab <- switch(type,hr="Hazard ratio",hazard="Hazard",surv="Survival",
                   sdiff="Survival difference",hdiff="Hazard difference",cumhaz="Cumulative hazard")
  xx <- attr(y,"newdata")
  xx <- eval(x@timeExpr,xx) # xx[,ncol(xx)]
  if (!add) matplot(xx, y, type="n", xlab=xlab, ylab=ylab, ...)
  if (ci) polygon(c(xx,rev(xx)), c(y[,2],rev(y[,3])), col=ci.col, border=ci.col)
  lines(xx,y[,1],col=line.col,lty=lty)
  if (rug) {
      Y <- x@y
      eventTimes <- Y[Y[,ncol(Y)]==1,ncol(Y)-1]
      rug(eventTimes,col=line.col)
    }
  return(invisible(y))
})

## penalised stpm2
setOldClass("gam")
setClass("pstpm2", representation(xlevels="list",
                                  contrasts="listOrNULL",
                                  terms="terms",
                                  logli="function",
                                  gam="gam",
                                  timeVar="character",
                                  timeExpr="nameOrcall",
                                  model.frame="list",
                                  call.formula="formula",
                                  x="matrix",
                                  xd="matrix",
                                  termsd="terms",
                                  Call="call",
                                  y="Surv"
                                  ),
         contains="mle2")
pstpm2 <- function(formula, data,
                     logH.args = NULL, logH.formula = NULL,
                     tvc = NULL, tvc.formula = NULL,
                     control = list(parscale = 0.1, maxit = 300), init = FALSE,
                     coxph.strata = NULL, weights = NULL, robust = FALSE, baseoff = FALSE,
                     bhazard = NULL, timeVar = NULL, sp=NULL, use.gr = TRUE,
                     contrasts = NULL, subset = NULL, ...)
  {
    ## set up the data
    ## ensure that data is a data frame
    data <- get_all_vars(formula, data)
    ## restrict to non-missing data (assumes na.action=na.omit)
    .include <- Reduce(`&`,
                       lapply(model.frame(formula, data, na.action=na.pass),
                              Negate(is.na)),
                       TRUE)
    data <- data[.include, , drop=FALSE]
    ##
    ## parse the function call
    Call <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "contrasts", "weights"),
               names(mf), 0L)
    mf <- mf[c(1L, m)]
    ##
    ## parse the event expression
    stopifnot(length(lhs(formula))>=2)
    eventExpr <- lhs(formula)[[length(lhs(formula))]]
    delayed <- length(lhs(formula))==4
    timeExpr <- lhs(formula)[[if (delayed) 3 else 2]] # expression
    if (is.null(timeVar))
      timeVar <- all.vars(timeExpr)
    time <- eval(timeExpr, data)
    if (delayed) {
      time0Expr <- lhs(formula)[[2]]
      time0 <- eval(time0Expr, data)
    }
    event <- eval(eventExpr,data)
    event <- event > min(event)
    ##
    ## set up the formulae
    if (is.null(logH.formula) && is.null(logH.args)) {
      logH.args$k <- -1
    }
    if (is.null(logH.formula))
      logH.formula <- as.formula(call("~",as.call(c(quote(s),call("log",timeExpr),
                                                    vector2call(logH.args)))))
    if (is.null(tvc.formula) && !is.null(tvc)) {
      tvc.formulas <-
        lapply(names(tvc), function(name)
               call(":",
                    as.name(name),
                    as.call(c(quote(s),
                              call("log",timeExpr),
                              vector2call(list(df=tvc[[name]]))))))
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
    ##
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
    data$logHhat <- pmax(-18,log(-log(Shat(coxph.obj))))
    ##
    ## initial values and object for lpmatrix predictions
    gam.call <- mf
    gam.call[[1L]] <- as.name("gam")
    gam.formula <- full.formula
    lhs(gam.formula) <- quote(logHhat) # new response
    gam.call$formula <- gam.formula
    gam.call$sp <- sp
    dataEvents <- data[event,]
    gam.call$data <- quote(dataEvents) # events only
    gam.obj <- eval(gam.call)
    if (!init) {
      init <- coef(gam.obj)
    }
    ##
    ## set up X, mf and wt
    X <- predict(gam.obj,data,type="lpmatrix")
    mt <- terms(gam.obj)
    mf <- model.frame(gam.obj)
    wt <- model.weights(gam.obj$model)
    if (is.null(wt)) wt <- rep(1,nrow(X))
    ##
    ## XD matrix
    lpfunc <- function(delta,fit,data,var) {
      data[[var]] <- data[[var]]+delta
      predict(fit,data,type="lpmatrix")
    }
    XD <- grad(lpfunc,0,gam.obj,data,timeVar)
    XD <- matrix(XD,nrow=nrow(X))
    ##
    bhazard <- substitute(bhazard)
    bhazard <- if (is.null(bhazard)) 0 else eval(bhazard,data,parent.frame())
    ## smoothing parameters
    if (is.null(sp))
      sp <- if(is.null(gam.obj$full.sp)) gam.obj$sp else gam.obj$full.sp
    ## penalty function
    pfun <- function(beta)
      sum(sapply(1:length(gam.obj$smooth),
                 function(i) {
                   smoother <- gam.obj$smooth[[i]]
                   betai <- beta[smoother$first.para:smoother$last.para]
                   sp[i]/2 * betai %*% smoother$S[[1]] %*% betai
                 }))
    dpfun <- function(beta) {
      deriv <- beta*0
      for (i in 1:length(gam.obj$smooth))
        {
          smoother <- gam.obj$smooth[[i]]
          ind <- smoother$first.para:smoother$last.para
          deriv[ind] <- sp[i] * smoother$S[[1]] %*% beta[ind]
        }
      return(deriv)
    }
    if (delayed && any(time0>0)) {
      ind <- time0>0
      data2 <- data[ind,,drop=FALSE] # data for delayed entry times
      X2 <- predict(gam.obj, data2, type="lpmatrix")
      wt2 <- wt[ind]
      gradnegll <- function(beta) {
        eta <- as.vector(X %*% beta)
        eta2 <- as.vector(X2 %*% beta)
        etaD <- as.vector(XD %*% beta)
        h <- etaD*exp(eta) + bhazard
        h[h<0] <- 1e-100
        g <- colSums(exp(eta)*wt*(-X + ifelse(event,1/h,0)*(XD + X*etaD)))-
          colSums(exp(eta2)*wt2*X2) - dpfun(beta)
        return(-g)
      }
      negll <- function(beta) {
        eta <- X %*% beta
        eta2 <- X2 %*% beta
        h <- (XD %*% beta)*exp(eta) + bhazard
        h[h<0] <- 1e-100
        ll <- sum(wt[event]*log(h[event])) +  sum(wt2*exp(eta2)) -
          sum(wt*exp(eta)) - pfun(beta)
        return(-ll)
      }
      logli <- function(beta) {
        eta <- X %*% beta
        eta2 <- X2 %*% beta
        h <- (XD %*% beta)*exp(eta) + bhazard
        h[h<0] <- 1e-100
        out <- exp(eta2) - exp(eta)
        out[event] <- out[event]+log(h[event])
        out <- out*wt
        return(out)
      }
    }
    else { # right censoring only
      gradnegll <- function(beta) {
        eta <- as.vector(X %*% beta)
        etaD <- as.vector(XD %*% beta)
        h <- etaD*exp(eta) + bhazard
        h[h<0] <- 1e-100
        g <- colSums(exp(eta)*wt*(-X + ifelse(event,1/h,0)*(XD + X*etaD))) - dpfun(beta)
        return(-g)
      }
      negll <- function(beta) {
        eta <- X %*% beta
        h <- (XD %*% beta)*exp(eta) + bhazard
        h[h<0] <- 1e-100
        ll <- sum(wt[event]*log(h[event])) - sum(wt*exp(eta)) - pfun(beta)
        return(-ll)
      }
      logli <- function(beta) {
        eta <- X %*% beta
        h <- (XD %*% beta)*exp(eta) + bhazard
        h[h<0] <- 1e-100
        out <- -exp(eta)
        out[event] <- out[event]+log(h[event])
        out <- out*wt
        return(out)
      }
    }
    ## MLE
    if (!is.null(control) && "parscale" %in% names(control)) {
      if (length(control$parscale)==1)
        control$parscale <- rep(control$parscale,length(init))
      if (is.null(names(control$parscale)))
        names(control$parscale) <- names(init)
    }
    parnames(negll) <- parnames(gradnegll) <- names(init)
    mle2 <- if (use.gr)
      mle2(negll,init,vecpar=TRUE, control=control, gr=gradnegll, ...)
    else mle2(negll,init,vecpar=TRUE, control=control, ...)
    ## mle2 <- mle2(negll,init,vecpar=TRUE, control=control, ...)
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
               data = data,
               formula = mle2@formula,
               optimizer = "optim",
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
               timeExpr = timeExpr,
               call.formula = formula,
               x = X,
               xd = XD,
               termsd = mt, # wrong!
               y = y)
    if (robust) # kludge
      out@vcov <- sandwich.stpm2(out)
    return(out)
  }
setMethod("predictnl", "pstpm2",
          function(object,fun,newdata=NULL,link=c("I","log","cloglog"),...)
  {
    link <- match.arg(link)
    invlinkf <- switch(link,I=I,log=exp,cloglog=cexpexp)
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
                   type=c("surv","cumhaz","hazard","hr","sdiff","hdiff","loghazard","link"),
                   grid=FALSE,seqLength=300,
                   se.fit=FALSE,link=NULL,exposed=incrVar(var),var,...)
  {
    ## exposed is a function that takes newdata and returns the revised newdata
    ## var is a string for a variable that defines a unit change in exposure
    local <-  function (object, newdata=NULL, type="surv", exposed)
      {
        tt <- object@terms 
        if (is.null(newdata)) {
          ##mm <- X <- model.matrix(object) # fails (missing timevar)
          X <- object@x
          XD <- object@xd
          ##y <- model.response(object@model.frame)
          y <- object@y
          time <- as.vector(y[,ncol(y)-1])
        }
        else {
          lpfunc <- function(delta,fit,data,var) {
            data[[var]] <- data[[var]]+delta
            predict(fit,data,type="lpmatrix")
          }
          X <- predict(object@gam, newdata, type="lpmatrix")
          XD <- grad(lpfunc,0,object@gam,newdata,object@timeVar)
          XD <- matrix(XD,nrow=nrow(X))
          ## resp <- attr(Terms, "variables")[attr(Terms, "response")] 
          ## similarly for the derivatives
          if (type %in% c("hazard","hr","sdiff","hdiff","loghazard")) {
            ## how to elegantly extract the time variable?
            timeExpr <- 
              lhs(object@call.formula)[[length(lhs(object@call.formula))-1]]
            time <- eval(timeExpr,newdata)
            ##
          }
          if (type %in% c("hr","sdiff","hdiff")) {
            if (missing(exposed))
              stop("exposed needs to be specified for type in ('hr','sdiff','hdiff')")
            newdata2 <- exposed(newdata)
            X2 <- predict(object@gam, newdata2, type="lpmatrix")
            XD2 <- grad(lpfunc,0,object@gam,newdata2,object@timeVar)
            XD2 <- matrix(XD,nrow=nrow(X))
          }
        }
        beta <- coef(object)
        cumHaz = exp(X %*% beta)
        Sigma = vcov(object)
        if (type=="link") { # delayed entry?
          return(X %*% beta)
        }
        if (type=="cumhaz") { # delayed entry?
          return(cumHaz)
        }
        if (type=="surv") { # delayed entry?
          return(exp(-cumHaz))
        }
        if (type=="sdiff")
          return(exp(-exp(X2 %*% beta)) - exp(-cumHaz))
        if (type=="hazard") {
          return((XD %*% beta)*cumHaz)
        }
        if (type=="loghazard") {
          return(log(XD %*% beta)+log(cumHaz))
        }
        if (type=="hdiff") {
          return((XD2 %*% beta)*exp(X2 %*% beta) - (XD %*% betaXD)/time*cumHaz)
        }
        if (type=="hr") {
          cumHazRatio = exp((X2 - X) %*% beta)
          return((XD2 %*% beta)/(XD %*% beta)*cumHazRatio)
        }
      }
    ##debug(local)
    type <- match.arg(type)
    if (is.null(newdata) && type %in% c("hr","sdiff","hdiff"))
      stop("Prediction using type in ('hr','sdiff','hdiff') requires newdata to be specified.")
    if (grid) {
      Y <- object@y
      event <- Y[,ncol(Y)]==1
      time <- object@data[[object@timeVar]]
      eventTimes <- time[event]
      X <- seq(min(eventTimes),max(eventTimes),length=seqLength)[-1]
      data.x <- data.frame(X)
      names(data.x) <- object@timeVar
      newdata <- merge(newdata,data.x)
    }
    pred <- if (!se.fit) {
      local(object,newdata,type=type,exposed=exposed,
            ...)
    }
    else {
      if (is.null(link))
        link <- switch(type,surv="cloglog",cumhaz="log",hazard="log",hr="log",sdiff="I",
                       hdiff="I",loghazard="I",link="I")
      predictnl(object,local,link=link,newdata=newdata,type=type,
                exposed=exposed,...) 
    }
    attr(pred,"newdata") <- newdata
    ##if (grid) cbind(newdata,as.data.frame(pred)) else pred
    return(pred)
  })
##`%c%` <- function(f,g) function(...) g(f(...)) # function composition
## to do:
## (*) Stata-compatible knots
setMethod("plot", signature(x="pstpm2", y="missing"),
          function(x,y,newdata,type="surv",
                      xlab=NULL,ylab=NULL,line.col=1,ci.col="grey",lty=par("lty"),
                      add=FALSE,ci=TRUE,rug=TRUE,
                      var=NULL,...) {
  y <- predict(x,newdata,type=type,var=var,grid=TRUE,se.fit=TRUE)
  if (is.null(xlab)) xlab <- deparse(x@timeExpr)
  if (is.null(ylab))
    ylab <- switch(type,hr="Hazard ratio",hazard="Hazard",surv="Survival",
                   sdiff="Survival difference",hdiff="Hazard difference",cumhaz="Cumulative hazard")
  xx <- attr(y,"newdata")
  xx <- eval(x@timeExpr,xx) # xx[,ncol(xx)]
  if (!add) matplot(xx, y, type="n", xlab=xlab, ylab=ylab, ...)
  if (ci) polygon(c(xx,rev(xx)), c(y[,2],rev(y[,3])), col=ci.col, border=ci.col)
  lines(xx,y[,1],col=line.col,lty=lty)
  if (rug) {
      Y <- x@y
      eventTimes <- Y[Y[,ncol(Y)]==1,ncol(Y)-1]
      rug(eventTimes,col=line.col)
    }
  return(invisible(y))
})

if (FALSE) {
try(suppressWarnings(detach("package:rstpm2",unload=TRUE)),silent=TRUE)
require(rstpm2)
data(brcancer)
system.time(fit2 <- stpm2(Surv(rectime/365,censrec==1)~hormon,data=brcancer,df=5))
system.time(fit3 <- pstpm2(Surv(rectime/365,censrec==1)~hormon,data=brcancer,use.gr=F))
plot(fit3,newdata=data.frame(hormon=0),type="hazard")
plot(fit2,newdata=data.frame(hormon=0),type="hazard",add=TRUE,ci=FALSE,rug=FALSE,
     line.col=2)
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
}


if (FALSE) { ##### examples #####
require(foreign)
if (FALSE) { # testing in open code
  install.packages("bbmle", repos="http://R-Forge.R-project.org")
  require(bbmle)
  brcancer=read.dta("brcancer.dta")
  brcancer=transform(brcancer,rate0=10^(-5+x1/100))
}
try(suppressWarnings(detach("package:bbmle",unload=TRUE)),silent=TRUE)

try(suppressWarnings(detach("package:rstpm2",unload=TRUE)),silent=TRUE)
require(rstpm2)
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

require(rstpm2)
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
require(lattice)
xyplot(Estimate ~ rectime, data=pred.all, group=hormon,type="l",xlab="Time")


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
brcancer10 = do.call("rbind",lapply(1:100,function(i) brcancer))
system.time(summary(fit <- stpm2(Surv(rectime,censrec==1)~hormon,df=3,data=brcancer10)))


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

} ## end of examples ##




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
