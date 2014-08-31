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
numDeltaMethod <- function(object,fun,...) {
  coef <- coef(object)
  est <- fun(coef,...)
  Sigma <- vcov(object)
  gd <- grad(fun,coef,...)
  ## se.est <- as.vector(sqrt(diag(t(gd) %*% Sigma %*% gd)))
  se.est <- as.vector(sqrt(colSums(gd* (Sigma %*% gd))))
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
  rr <- t(grad(obj@logli,coef(obj)))
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
setClassUnion("nameOrcallOrNULL",c("name","call","NULL"))
##setClassUnion("numericOrNULL",c("numeric","NULL"))
setOldClass("Surv")
setOldClass("lm")

## re-factored stpm2
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
                                 model.frame="list",
                                 call.formula="formula",
                                 x="matrix",
                                 xd="matrix",
                                 termsd="terms",
                                 Call="call",
                                 y="Surv"
                                 ),
         contains="mle2")
## check: weights
stpm2 <- function(formula, data,
                     df = 3, cure = FALSE, logH.args = NULL, logH.formula = NULL,
                     tvc = NULL, tvc.formula = NULL,
                     control = list(parscale = 0.1, maxit = 300), init = FALSE,
                     coxph.strata = NULL, weights = NULL, robust = FALSE, baseoff = FALSE,
                     bhazard = NULL, timeVar = "", time0Var = "", use.gr = TRUE, use.rcpp= TRUE,
                  reltol=1.0e-8, 
                     contrasts = NULL, subset = NULL, ...)
  {
    ## parse the event expression
    stopifnot(length(lhs(formula))>=2)
    eventExpr <- lhs(formula)[[length(lhs(formula))]]
    delayed <- length(lhs(formula))==4
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
    ##
    ## set up the data
    ## ensure that data is a data frame
    data <- get_all_vars(full.formula, data)
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
    ## get variables
    time <- eval(timeExpr, data)
    time0Expr <- NULL # initialise
    if (delayed) {
      time0Expr <- lhs(formula)[[2]]
      if (time0Var == "")
        time0Var <- all.vars(time0Expr)
      time0 <- eval(time0Expr, data)
    }
    event <- eval(eventExpr,data)
    event <- event > min(event)
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
    if (is.logical(init) && !init) {
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
    lpfunc <- function(delta,fit,dataset,var) {
      dataset[[var]] <- dataset[[var]]+delta
      lpmatrix.lm(fit,dataset)
    }
    XD <- grad(lpfunc,0,lm.obj,data,timeVar)
    XD <- matrix(XD,nrow=nrow(X))
    ##
    bhazard <- substitute(bhazard)
    bhazard <- if (is.null(bhazard)) rep(0,nrow(X)) else eval(bhazard,data,parent.frame())
    if (delayed && all(time0==0)) delayed <- FALSE
    if (delayed) {
        ind <- time0>0
        data0 <- data[ind,,drop=FALSE] # data for delayed entry times
        X0 <- lpmatrix.lm(lm.obj, data0)
        wt0 <- wt[ind]
    } else {
        X0 <- wt0 <- NULL
    }
    gradnegll <- function(beta) {
        eta <- as.vector(X %*% beta)
        etaD <- as.vector(XD %*% beta)
        h <- etaD*exp(eta) + bhazard
        ## h[h<0] <- 1e-100
        g <- colSums(exp(eta)*wt*(-X + ifelse(event,1/h,0)*(XD + X*etaD)))
        if (delayed) {
            eta0 <- as.vector(X0 %*% beta)
            g <- g + colSums(exp(eta0)*wt0*X0)
        }
        return(-g)
      }
      negll <- function(beta) {
        eta <- X %*% beta
        h <- (XD %*% beta)*exp(eta) + bhazard
        ## h[h<0] <- 1e-100
        ll <- sum(wt[event]*log(h[event])) - sum(wt*exp(eta))
        if (delayed) {
            eta0 <- X0 %*% beta
            ll <- ll + sum(wt0*exp(eta0))
        }
        return(-ll)
      }
    logli <- function(beta) {
        eta <- X %*% beta
        h <- (XD %*% beta)*exp(eta) + bhazard
        ##  h[h<0] <- 1e-100
        out <- - exp(eta)
        out[event] <- out[event]+log(h[event])
        if (delayed) {
            eta0 <- X0 %*% beta
            out[ind] <- out[ind] + exp(eta0)
        }
        out <- out*wt
        return(out)
    }
    if (!is.null(control) && "parscale" %in% names(control)) {
      if (length(control$parscale)==1)
        control$parscale <- rep(control$parscale,length(init))
      if (is.null(names(control$parscale)))
        names(control$parscale) <- names(init)
    }
    parnames(negll) <- parnames(gradnegll) <- names(init)
    rcpp_stpm2 <- function() {
        stopifnot(!delayed)
        parscale <- if (!is.null(control$parscale)) control$parscale else rep(1,length(init))
        names(parscale) <- names(init)
        .Call("optim_stpm2",list(init=init,X=X,XD=XD,bhazard=bhazard,wt=wt,event=ifelse(event,1,0),
              delayed=if (delayed) 1 else 0, X0=X0, wt0=wt0, parscale=parscale, reltol=reltol),
              package="rstpm2")
    }
    analyticalHessian <- function(beta) {
        mult <- function(X1,X2,w) {
            out <- matrix(0,ncol(X1),ncol(X1))
            for (k in (1:nrow(X1))[w != 0])
                out <- out + (X1[k,] %*% t(X2[k,])) * w[k]
            out
        }
        eta <- X %*% beta
        etaD <- XD %*% beta
        expEta <- exp(eta)
        h <- etaD*expEta + bhazard
        ## case: bhazard=0
        ## hess <- - mult(X,X,expEta*wt) - mult(XD,XD,event*wt/(etaD^2))
        w1 <- event/h*wt*expEta
        w2 <- event/h^2*wt*expEta^2
        hess <- - mult(X,X,expEta*wt) +
            mult(X,X,w1*etaD) + mult(X,XD,w1) + mult(XD,X,w1) -
                (mult(X,X,w2*etaD^2) + mult(XD,X,w2*etaD) + mult(X,XD,w2*etaD) + mult(XD,XD,w2))
        ## print(solve(-hess))
        hess
    }
    ## MLE
    if (use.rcpp) {
        fit <- rcpp_stpm2()
        coef <- as.vector(fit$coef)
        hessian <- fit$hessian
        names(coef) <- rownames(hessian) <- colnames(hessian) <- names(init)
        mle2 <- mle2(negll, coef, vecpar=TRUE, control=control, gr=gradnegll, ..., eval.only=TRUE)
        mle2@vcov <- solve(hessian)
        mle2@details$convergence <- fit$fail
    } else {
        mle2 <- if (use.gr)
            mle2(negll,init,vecpar=TRUE, control=control, gr=gradnegll, ...)
        else mle2(negll,init,vecpar=TRUE, control=control, ...)
    }
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
               time0Var = time0Var,
               timeExpr = timeExpr,
               time0Expr = time0Expr,
               delayed = delayed,
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
            ## timeExpr <- 
            ##   lhs(object@call.formula)[[length(lhs(object@call.formula))-1]]
            time <- eval(object@timeExpr,newdata)
            ##
          }
          if (object@delayed) {
            newdata0 <- newdata
            newdata0[[object@timeVar]] <- newdata[[object@time0Var]]
            X0 <- lpmatrix.lm(object@lm, newdata0)
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
          if (object@delayed) return(cumHaz - (X0 %*% beta)) else return(cumHaz)
          ##return(cumHaz)
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
          return((XD2 %*% beta)*exp(X2 %*% beta) - (XD %*% beta)/time*cumHaz)
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
                      add=FALSE,ci=!add,rug=!add,
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
                                  timeExpr="nameOrcall",
                                  like="function",
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
                   coxph.strata = NULL, nStrata=5, weights = NULL, robust = FALSE, baseoff = FALSE,
                   bhazard = NULL, timeVar = NULL, sp=NULL, use.gr = TRUE, use.rcpp = TRUE, criterion=c("GCV","BIC"), penalty = c("logH","h"), smoother.parameters = NULL,
                   alpha=switch(criterion,BIC=1,GCV=1.4), sp.init=NULL, trace = 0,
                   reltol = list(search = 1.0e-6, final = 1.0e-8),
                   contrasts = NULL, subset = NULL, ...)
  {
    ## set up the data
    ## ensure that data is a data frame
      temp.formula <- formula
      if (!is.null(logH.formula)) rhs(temp.formula) <-rhs(temp.formula) %call+% rhs(logH.formula)
      if (!is.null(tvc.formula)) rhs(temp.formula) <-rhs(temp.formula) %call+% rhs(tvc.formula)
      raw.data <- data
    data <- get_all_vars(temp.formula, raw.data)
    criterion <- match.arg(criterion)
    penalty <- match.arg(penalty)
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
                              vector2call(list(k=tvc[[name]]))))))
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
    coxph.call$data <- quote(coxph.data)
    coxph.data <- data
    if (!is.null(coxph.strata)) {
        coxph.formula <- coxph.call$formula
        rhs(coxph.formula) <- rhs(formula) %call+% call("strata",coxph.strata)
        coxph.call$formula <- coxph.formula
    }
    if (nStrata) {
        coxph.strata <- substitute(GrOuP)
        coxph.data[["GrOuP"]] <- sample(1:nrow(data) %% nStrata ,nrow(data))
        coxph.formula <- coxph.call$formula
        rhs(coxph.formula) <- rhs(formula) %call+% call("strata",coxph.strata)
        coxph.call$formula <- coxph.formula
    }
    coxph.call$model <- TRUE
    ## coxph.obj <- eval(coxph.call, envir=parent.frame())
    coxph.obj <- eval(coxph.call, coxph.data)
    y <- model.extract(model.frame(coxph.obj),"response")
    data$logHhat <- pmax(-18,log(-log(Shat(coxph.obj))))
    ##
    ## initial values and object for lpmatrix predictions
    gam.call <- mf
    gam.call[[1L]] <- as.name("gam")
    gam.formula <- full.formula
    lhs(gam.formula) <- quote(logHhat) # new response
    gam.call$formula <- gam.formula
    ## gam.call$sp <- if (is.list(sp)) sp[[floor(length(sp)/2)]] else sp
    gam.call$sp <- sp
    if (is.null(sp) && !is.null(sp.init))
        gam.call$sp <- sp.init
    dataEvents <- data[event,]
    gam.call$data <- quote(dataEvents) # events only
    gam.obj <- eval(gam.call)
    if (is.logical(init) && !init) {
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
    ## lpfunc <- function(delta,fit,data,var) {
    ##   data[[var]] <- data[[var]]+delta
    ##   predict(fit,data,type="lpmatrix")
    ## }
    ## XD2 <- grad(lpfunc,0,gam.obj,data,timeVar)
    ## XD2 <- matrix(XD,nrow=nrow(X))
    lpfunc <- function(x,...) {
      newdata <- data
      newdata[[timeVar]] <- x
      predict(gam.obj,newdata,type="lpmatrix")
    }
    XD <- grad1(lpfunc,data[[timeVar]])    
    ##
    bhazard <- substitute(bhazard)
    bhazard <- if (is.null(bhazard)) rep(0,nrow(X)) else eval(bhazard,data,parent.frame())
    ## smoothing parameters
    if (no.sp <- is.null(sp)) {
        sp <- if(is.null(gam.obj$full.sp)) gam.obj$sp else gam.obj$full.sp
        if (!is.null(sp.init)) sp <- sp.init
    }
    ## penalty function
    pfun <- function(beta,sp)
      sum(sapply(1:length(gam.obj$smooth),
                 function(i) {
                   smoother <- gam.obj$smooth[[i]]
                   betai <- beta[smoother$first.para:smoother$last.para]
                   sp[i]/2 * betai %*% smoother$S[[1]] %*% betai
                 }))
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
    if (delayed && all(time0==0)) delayed <- FALSE
    if (delayed) {
        ind <- time0>0
        data0 <- data[ind,,drop=FALSE] # data for delayed entry times
        X0 <- predict(gam.obj, data0, type="lpmatrix")
        wt0 <- wt[ind]
    } else {
        X0 <- wt0 <- NULL
    }
    gradnegllsp <- function(beta,sp) {
        eta <- as.vector(X %*% beta)
        etaD <- as.vector(XD %*% beta)
        h <- etaD*exp(eta) + bhazard
        g <- colSums(exp(eta)*wt*(-X + ifelse(event,1/h,0)*(XD + X*etaD))) - dpfun(beta,sp)
        if (delayed) {
            eta0 <- as.vector(X0 %*% beta)
            g <- g + colSums(exp(eta0)*wt0*X0)
        }
        return(-g)
      }
    negllsp <- function(beta,sp) {
        eta <- X %*% beta
        h <- (XD %*% beta)*exp(eta) + bhazard
        ll <- sum(wt[event]*log(h[event])) - sum(wt*exp(eta)) - pfun(beta,sp)
        if (delayed) {
            eta0 <- X0 %*% beta
            ll <- ll + sum(wt0*exp(eta0))
        }
        return(-ll)
      }
    logli <- function(beta) {
        eta <- X %*% beta
        h <- (XD %*% beta)*exp(eta) + bhazard
        out <- - exp(eta)
        out[event] <- out[event]+log(h[event])
        if (delayed) {
            eta0 <- X0 %*% beta
            out[ind] <- out[ind] + exp(eta0)
        }
        out <- out*wt
        return(out)
    }
    like <- function(beta) {
        eta <- X %*% beta
        h <- (XD %*% beta)*exp(eta) + bhazard
        ll <- sum(wt[event]*log(h[event])) - sum(wt*exp(eta))
        if (delayed) {
            eta0 <- X0 %*% beta
            ll <- ll + sum(wt0*exp(eta0))
        }
        return(ll)
      }
    ## MLE
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
    rcpp_optim <- function() {
        ## stopifnot(!delayed)
        args <- list(init=init,X=X,XD=XD,bhazard=bhazard,wt=wt,event=ifelse(event,1,0),
                     delayed=if (delayed) 1 else 0, X0=X0, wt0=wt0, parscale=control$parscale,
                     smooth=if(penalty == "logH") gam.obj$smooth else design,
                     sp=sp, reltol_search=reltol$search, reltol=reltol$final, trace=trace,
                     alpha=alpha,criterion=switch(criterion,GCV=1,BIC=2))
        if (!no.sp) { # fixed sp as specified
          if (penalty == "logH")
              .Call("optim_pstpm2LogH_fixedsp", args, package = "rstpm2") else
          .Call("optim_pstpm2Haz_fixedsp", args, package = "rstpm2")
        }
        else if (length(sp)>1) {
            if (penalty == "logH")
            .Call("optim_pstpm2LogH_multivariate", args, package = "rstpm2") else
            .Call("optim_pstpm2Haz_multivariate", args, package = "rstpm2")
        } else {
            if (penalty == "logH")
            .Call("optim_pstpm2LogH_first", args, package = "rstpm2") else
            .Call("optim_pstpm2Haz_first", args, package = "rstpm2")
        }
    }
    if (use.rcpp) {
        fit <- rcpp_optim()
        fit$coef <- as.vector(fit$coef)
        names(fit$coef) <- names(init)
        init <- fit$coef
        if (!no.sp) sp <- fit$sp # ignores alpha
    }
    negll <- function(beta) negllsp(beta,sp)
    gradnegll <- function(beta) gradnegllsp(beta,sp)
    parnames(negll) <- parnames(gradnegll) <- names(init)
    if (use.rcpp) {
        mle2 <- if (use.gr) {
            mle2(negll,init,vecpar=TRUE, control=control, gr=gradnegll, eval.only=TRUE, ...)
        } else mle2(negll,init,vecpar=TRUE, control=control, eval.only=TRUE, ...)
        mle2@details$hessian <- fit$hessian
        mle2@vcov <- solve(fit$hessian)
        mle2@details$convergence <- 0
    } else {
        mle2 <- if (use.gr) {
            mle2(negll,init,vecpar=TRUE, control=control, gr=gradnegll, ...)
        } else mle2(negll,init,vecpar=TRUE, control=control, ...)
        if (any(is.na(mle2@vcov)))
            mle2@vcov <- solve(optimHess(coef(mle2),negll,gradnegll))
    }
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
                   timeExpr = timeExpr,
                   like = like,
                   call.formula = formula,
                   x = X,
                   xd = XD,
                   termsd = mt, # wrong!
                   y = y)
    if (robust) # kludge
        out@vcov <- sandwich.stpm2(out)
    return(out)
  }

########GCV##############
##require(numDeriv)
## now fit a penalised stpm2 model
##pstpm2.fit <- pstpm2(formula,data)
## log likelihood and penalized log likelihood

##GCV###
gcv<-function(pstpm2.fit){
  like<-pstpm2.fit@like
  Hl<-numDeriv::hessian(like,coef(pstpm2.fit))
  if (any(is.na(Hl)))
      Hl <- optimHess(coef(pstpm2.fit), like)
  Hinv<- -vcov(pstpm2.fit)
  trace<-sum(diag(Hinv%*%Hl))
  l<-like(coef(pstpm2.fit))
  structure(-l+trace,negll=-l,trace=trace)
}
###AICC
aicc<-function(pstpm2.fit,nn){
  like<-pstpm2.fit@like
  Hl<-numDeriv::hessian(like,coef(pstpm2.fit))
  Hinv<--vcov(pstpm2.fit)
  trace<-sum(diag(Hinv%*%Hl))
  ll<-like(coef(pstpm2.fit))
return(-2*ll+2*trace*nn/(nn-trace-1))
}
###BIC
bic<-function(pstpm2.fit,nn){
  like<-pstpm2.fit@like
  Hl<-numDeriv::hessian(like,coef(pstpm2.fit))
  Hinv<--vcov(pstpm2.fit)
 trace<-sum(diag(Hinv%*% Hl))
ll<-like(coef(pstpm2.fit))
return(-2*ll+trace*log(nn))
}
###GCVC
gcvc<-function(pstpm2.fit,nn){
  like<-pstpm2.fit@like
  Hl<-numDeriv::hessian(like,coef(pstpm2.fit))
  Hinv<-vcov(pstpm2.fit)
  trace<-sum(diag(Hinv %*% Hl))
  ll<-like(coef(pstpm2.fit))
return(-2*ll-2*nn*log(1-trace/nn))
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
          X <- predict(object@gam, newdata, type="lpmatrix")
          ## lpfunc <- function(delta,fit,data,var) {
          ##   data[[var]] <- data[[var]]+delta
          ##   predict(fit,data,type="lpmatrix")
          ## }
          ## XD <- grad(lpfunc,0,object@gam,newdata,object@timeVar)
          ## XD <- matrix(XD,nrow=nrow(X))
          lpfunc <- function(x,...) {
            newdata2 <- newdata
            newdata2[[object@timeVar]] <- x
            predict(object@gam,newdata2,type="lpmatrix")
          }
          XD <- grad1(lpfunc,newdata[[object@timeVar]])    
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
                   lwd=par("lwd"),
                   add=FALSE,ci=!add,rug=!add,
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
  lines(xx,y[,1],col=line.col,lty=lty,lwd=lwd)
  if (rug) {
      Y <- x@y
      eventTimes <- Y[Y[,ncol(Y)]==1,ncol(Y)-1]
      rug(eventTimes,col=line.col)
    }
  return(invisible(y))
})
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
## if (FALSE) {
## #####load data####
## load("brcancer.rda")
## data(brcancer)
## brcancer$recyear <- brcancer$rectime/365
## ####model fit###
## opt.fit(Surv(recyear,censrec==1)~hormon,data=brcancer,
##         logH.formula=~s(recyear), sp.low=10^-4,sp.upp=4000,
##         num.sp=30,timeVar = NULL)
## }
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
    bhazard <- if (is.null(bhazard)) rep(0,nrow(X)) else bhazard[event] # crude
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


