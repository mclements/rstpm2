

nsxD <- 
function (x, df = NULL, knots = NULL, intercept = FALSE, Boundary.knots = range(x), 
    derivs = if (cure) c(2, 1) else c(2, 2), log = FALSE, centre = FALSE, 
    cure = FALSE, stata.stpm2.compatible = FALSE) 
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
            else stats::quantile(x[!outside], round(knots, 2), 
                type = 2)
        }
    }
    else nIknots <- length(knots)
    Aknots <- sort(c(rep(Boundary.knots, 4L), knots))
    if (any(outside)) {
        basis <- array(0, c(length(x), nIknots + 4L))
        if (any(ol)) {
            k.pivot <- Boundary.knots[1L]
            tt <- spline.des(Aknots, rep(k.pivot, sum(ol)), 4, 1)$design
            basis[ol, ] <- tt
        }
        if (any(or)) {
            k.pivot <- Boundary.knots[2L]
            tt <- spline.des(Aknots, rep(k.pivot, sum(or)), 4, 1)$design
            basis[or, ] <- tt
        }
        if (any(inside <- !outside)) 
            basis[inside, ] <- spline.des(Aknots, x[inside], 
                4, 1)$design
    }
    else basis <- spline.des(Aknots, x, 4, 1)$design
    const <- splineDesign(Aknots, rep(Boundary.knots, 3 - derivs), 
        4, c(derivs[1]:2, derivs[2]:2))
    if (!intercept) {
        const <- const[, -1, drop = FALSE]
        basis <- basis[, -1, drop = FALSE]
    }
    qr.const <- qr(t(const))
    q.const <- qr.Q(qr.const, complete=TRUE)[, -(1L:2L), drop = FALSE] # NEW
    basis <- as.matrix((t(qr.qty(qr.const, t(basis))))[, -(1L:nrow(const)), 
        drop = FALSE])
    n.col <- ncol(basis)
    if (nas) {
        nmat <- matrix(NA, length(nax), n.col)
        nmat[!nax, ] <- basis
        basis <- nmat
    }
    dimnames(basis) <- list(nx, 1L:n.col)
    if (centre) {
        centreBasis <- nsx(centre, knots = if (is.null(knots)) 
            numeric(0)
        else knots, Boundary.knots = Boundary.knots, intercept = intercept, 
            derivs = derivs, centre = FALSE, log = log)
        oldAttributes <- attributes(basis)
        basis <- t(apply(basis, 1, function(x) x - centreBasis))
        attributes(basis) <- oldAttributes
    }
    a <- list(degree = 3, knots = if (is.null(knots)) numeric(0) else knots, 
        Boundary.knots = Boundary.knots, intercept = intercept, 
        derivs = derivs, centre = centre, log = log, q.const = q.const)
    attributes(basis) <- c(attributes(basis), a)
    class(basis) <- c("nsxD", "basis", "matrix")
    basis
}
makepredictcall.nsxD <- 
function (var, call) 
{
    if (as.character(call)[1L] != "nsxD") 
        return(call)
    at <- attributes(var)[c("knots", "Boundary.knots", "intercept",
                            "derivs", "centre", "log")]
    xxx <- call[1L:2]
    xxx[names(at)] <- at
    xxx
}
predict.nsxD <- 
function (object, newx, ...) 
{
    if (missing(newx)) 
        return(object)
    a <- c(list(x = newx), attributes(object)[c("knots", "Boundary.knots", 
        "intercept", "derivs", "centre", "log")])
    do.call("nsxD", a)
}


nsxDD <- 
function (x, df = NULL, knots = NULL, intercept = FALSE, Boundary.knots = range(x), 
    derivs = if (cure) c(2, 1) else c(2, 2), log = FALSE, centre = FALSE, 
    cure = FALSE, stata.stpm2.compatible = FALSE) 
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
            else stats::quantile(x[!outside], round(knots, 2), 
                type = 2)
        }
    }
    else nIknots <- length(knots)
    Aknots <- sort(c(rep(Boundary.knots, 4L), knots))
    if (any(outside)) {
        basis <- array(0, c(length(x), nIknots + 4L))
        if (any(ol)) {
            basis[ol, ] <- 0
        }
        if (any(or)) {
            basis[or, ] <- 0
        }
        if (any(inside <- !outside)) 
            basis[inside, ] <- spline.des(Aknots, x[inside], 
                4, 2)$design
    }
    else basis <- spline.des(Aknots, x, 4, 2)$design
    const <- splineDesign(Aknots, rep(Boundary.knots, 3 - derivs), 
        4, c(derivs[1]:2, derivs[2]:2))
    if (!intercept) {
        const <- const[, -1, drop = FALSE]
        basis <- basis[, -1, drop = FALSE]
    }
    qr.const <- qr(t(const))
    q.const <- qr.Q(qr.const, complete=TRUE)[, -(1L:2L), drop = FALSE] # NEW
    basis <- as.matrix((t(qr.qty(qr.const, t(basis))))[, -(1L:nrow(const)), 
        drop = FALSE])
    n.col <- ncol(basis)
    if (nas) {
        nmat <- matrix(NA, length(nax), n.col)
        nmat[!nax, ] <- basis
        basis <- nmat
    }
    dimnames(basis) <- list(nx, 1L:n.col)
    if (centre) {
        centreBasis <- nsx(centre, knots = if (is.null(knots)) 
            numeric(0)
        else knots, Boundary.knots = Boundary.knots, intercept = intercept, 
            derivs = derivs, centre = FALSE, log = log)
        oldAttributes <- attributes(basis)
        basis <- t(apply(basis, 1, function(x) x - centreBasis))
        attributes(basis) <- oldAttributes
    }
    a <- list(degree = 3, knots = if (is.null(knots)) numeric(0) else knots, 
        Boundary.knots = Boundary.knots, intercept = intercept, 
        derivs = derivs, centre = centre, log = log, q.const = q.const)
    attributes(basis) <- c(attributes(basis), a)
    class(basis) <- c("nsxDD", "basis", "matrix")
    basis
}
makepredictcall.nsxDD <- 
function (var, call) 
{
    if (as.character(call)[1L] != "nsxDD") 
        return(call)
    at <- attributes(var)[c("knots", "Boundary.knots", "intercept",
                            "derivs", "centre", "log")]
    xxx <- call[1L:2]
    xxx[names(at)] <- at
    xxx
}
predict.nsxDD <- 
function (object, newx, ...) 
{
    if (missing(newx)) 
        return(object)
    a <- c(list(x = newx), attributes(object)[c("knots", "Boundary.knots", 
        "intercept", "derivs", "centre", "log")])
    do.call("nsxDD", a)
}

## test nsxD and nsxDD
if (FALSE) {
    zeros <- function(mat,rows=1:nrow(mat),cols=1:ncol(mat)) "[<-"(mat,rows,cols,0)
    tm <- as.numeric(3:5)
    tm2 <- as.numeric(0:11)
    y <- rnorm(length(tm))
    lm1 <- lm(y~nsx(tm,df=4))
    lmD1 <- lm(y~nsxD(tm,df=4)-1)
    lmDD1 <- lm(y~nsxDD(tm,df=4)-1)
    eps <- 1e-5
    (lpmatrix.lm(lm1,newdata=data.frame(tm=tm2+eps)) - 
     lpmatrix.lm(lm1,newdata=data.frame(tm=tm2-eps)))/(2*eps) -
        cbind(0,lpmatrix.lm(lmD1,newdata=data.frame(tm=tm2))) # ok
    (lpmatrix.lm(lmD1,newdata=data.frame(tm=tm2+eps)) - 
     lpmatrix.lm(lmD1,newdata=data.frame(tm=tm2-eps)))/(2*eps) -
        lpmatrix.lm(lmDD1,newdata=data.frame(tm=tm2)) # ok
}

S0hat <- function(obj)
  {
    ## predicted survival for individuals (adjusted for covariates)
    newobj = survfit(obj,se.fit=FALSE)
    surv = newobj$surv
    surv[match(obj$y[,ncol(obj$y)-1],newobj$time)]
  }

## general link functions
setClass("aft", representation(args="list"), contains="mle2")

aft <- function(formula, data, smooth.formula = NULL, df = 3,
                 control = list(parscale = 1, maxit = 1000), init = NULL,
                 weights = NULL, 
                 timeVar = "", time0Var = "", log.time.transform=TRUE,
                 reltol=1.0e-8, trace = 0,
                 contrasts = NULL, subset = NULL, use.gr = TRUE, ...) {
    ## parse the event expression
    eventInstance <- eval(lhs(formula),envir=data)
    stopifnot(length(lhs(formula))>=2)
    eventExpr <- lhs(formula)[[length(lhs(formula))]]
    delayed <- length(lhs(formula))>=4 # indicator for multiple times (cf. strictly delayed)
    surv.type <- attr(eventInstance,"type")
    if (surv.type %in% c("interval","interval2","left","mstate"))
        stop("stpm2 not implemented for Surv type ",surv.type,".")
    counting <- attr(eventInstance,"type") == "counting"
    ## interval <- attr(eventInstance,"type") == "interval"
    timeExpr <- lhs(formula)[[if (delayed) 3 else 2]] # expression
    if (timeVar == "")
        timeVar <- all.vars(timeExpr)
    ## set up the formulae
    full.formula <- formula
    if (!is.null(smooth.formula))
        rhs(full.formula) <- rhs(formula) %call+% rhs(smooth.formula)
    rhs(full.formula) <- rhs(full.formula) %call+% quote(0)
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
    } else {
        time0 <- NULL
    }
    event <- eval(eventExpr,data)
    ## if all the events are the same, we assume that they are all events, else events are those greater than min(event)
    event <- if (length(unique(event))==1) rep(TRUE, length(event)) else event <- event > min(event)
    ## setup for initial values
    ## Cox regression
    coxph.call <- mf
    coxph.call[[1L]] <- as.name("coxph")
    coxph.call$model <- TRUE
    coxph.obj <- eval(coxph.call, envir=parent.frame())
    y <- model.extract(model.frame(coxph.obj),"response")
    data$logHhat <- pmax(-18,log(-log(S0hat(coxph.obj))))
    ##
    ## Weibull regression
    if (delayed) {
        if (requireNamespace("eha", quietly = TRUE)) {
            survreg1 <- eha::aftreg(formula, data)
            coef1 <- coef(survreg1)
            coef1 <- coef1[1:(length(coef1)-2)]
        } else coef1 <- rep(0,ncol(X))
    } else {
        survreg1 <- survival::survreg(formula, data)
        coef1 <- coef(survreg1)
        coef1 <- coef1[-1] # assumes intercept included in the formula; ignores smooth.formula
    }
    ## pred1 <- predict(survreg1)
    data$logtstar <- log(time)    
    ## data$logtstar <- log(time/pred1)    
    ## initial values and object for lpmatrix predictions
    lm.call <- mf
    lm.call[[1L]] <- as.name("lm")
    lm.formula <- full.formula
    lhs(lm.formula) <- quote(logtstar) # new response
    lm.call$formula <- lm.formula
    dataEvents <- data[event,]
    lm.call$data <- quote(dataEvents) # events only
    lm.obj <- eval(lm.call)
    coef1b <- coef(lm.obj)
    if (names(coef1b)[1]=="(Intercept)") coef1b <- coef1b[-1] # ???
    ## if (is.null(init)) {
    ##   init <- coef(lm.obj)
    ## }
    lm0.obj <- lm(logHhat~nsx(logtstar,df,intercept=TRUE)-1,dataEvents)
    ## lm0D.obj <- lm(logHhat~nsxD(logtstar,df,intercept=TRUE)-1,dataEvents)
    coef0 <- coef(lm0.obj) # log-log baseline
    ## design information for baseline survival
    design <- nsx(dataEvents$logtstar, df=df, intercept=TRUE)
    designD <- nsxD(dataEvents$logtstar, df=df, intercept=TRUE)
    designDD <- nsxDD(dataEvents$logtstar, df=df, intercept=TRUE)
    ##
    ## set up mf and wt
    mt <- terms(lm.obj)
    mf <- model.frame(lm.obj)
    ## wt <- model.weights(lm.obj$model)
    wt <- if (is.null(substitute(weights))) rep(1,nrow(data)) else eval(substitute(weights),data,parent.frame())
    ##
    ## XD matrix
    lpfunc <- function(x,fit,data,var) {
      data[[var]] <- x
      lpmatrix.lm(fit,data)
    }
    ##
    ## initialise values
    ind0 <- FALSE
    map0 <- 0L
    which0 <- 0
    wt0 <- 0
    ## ttype <- 0
    ## surv.type %in% c("right","counting")
    X <- lpmatrix.lm(lm.obj,data)
    XD <- grad1(lpfunc,data[[timeVar]],lm.obj,data,timeVar,log.transform=log.time.transform)
    XD <- matrix(XD,nrow=nrow(X))
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
        X0 <- lpmatrix.lm(lm.obj, data0)
        wt0 <- wt[ind0]
        rm(data0)
    }
    if (ncol(X)>length(coef1)) {
        coef1 <- c(coef1,rep(0,ncol(X) - length(coef1)))
        names(coef1) <- names(coef1b)
        }
    init <- c(coef1,coef0)
    if (!is.null(control) && "parscale" %in% names(control)) {
      if (length(control$parscale)==1)
        control$parscale <- rep(control$parscale,length(init))
      if (is.null(names(control$parscale)))
        names(control$parscale) <- names(init)
    }
    parscale <- rep(if (is.null(control$parscale)) 1 else control$parscale,length=length(init))
    names(parscale) <- names(init)
    args <- list(init=init,X=X,XD=XD,wt=wt,event=ifelse(event,1,0),time=time,y=y,
                 timeVar=timeVar,timeExpr=timeExpr,terms=mt,
                 delayed=delayed, X0=X0, wt0=wt0, parscale=parscale, reltol=reltol,
                 time0=if (delayed) time0[time0>0] else NULL, log.time.transform=log.time.transform,
                 trace = as.integer(trace), map0 = map0 - 1L, ind0 = ind0, which0 = which0 - 1L,
                 boundaryKnots=attr(design,"Boundary.knots"), q.const=t(attr(design,"q.const")),
                 interiorKnots=attr(design,"knots"), design=design, designD=designD,
                 designDD=designDD,
                 data=data, lm.obj = lm.obj, return_type="optim")
    negll <- function(beta) {
        localargs <- args
        localargs$return_type <- "objective"
        localargs$init <- beta
        return(.Call("aft_model_output", localargs, PACKAGE="rstpm2"))
    }
    gradient <- function(beta) {
        localargs <- args
        localargs$return_type <- "gradient"
        localargs$init <- beta
        return(.Call("aft_model_output", localargs, PACKAGE="rstpm2"))
    }
    negll.slow <- function(betafull) {
        beta <- betafull[1:ncol(args$X)]
        betas <- betafull[-(1:ncol(args$X))]
        time <- args$time
        eta <- as.vector(args$X %*% beta)
        etaD <- as.vector(args$XD %*% beta)
        logtstar <- log(args$time) - eta
        etas <- as.vector(predict(args$design,logtstar) %*% betas)
        etaDs <- as.vector(predict(args$designD,logtstar) %*% betas)
        ## fix bounds on etaDs
        eps <- etaDs*0. + 1e-8;
        pen <- sum(pmin(etaDs,eps)*pmin(etaDs,eps))
        etaDs <- pmax(etaDs, eps)
        ## fix bounds on etaD
        pen = pen + sum(pmin(1/time-etaD,eps)*pmin(1/time-etaD,eps))
        etaD <- 1/time - pmax(1/time-etaD, eps);
        logh <- etas + log(etaDs) + log(1/time -etaD)
        H <- exp(etas)
        pen - (sum(logh*event) - sum(H))
    }
    neglli <- function(betafull) {
        beta <- betafull[1:ncol(args$X)]
        betas <- betafull[-(1:ncol(args$X))]
        time <- args$time
        eta <- as.vector(args$X %*% beta)
        etaD <- as.vector(args$XD %*% beta)
        logtstar <- log(args$time) - eta
        etas <- as.vector(predict(args$design,logtstar) %*% betas)
        etaDs <- as.vector(predict(args$designD,logtstar) %*% betas)
        ## fix bounds on etaDs
        eps <- etaDs*0. + 1e-8;
        pen <- pmin(etaDs,eps)*pmin(etaDs,eps)
        etaDs <- pmax(etaDs, eps)
        ## fix bounds on etaD
        pen <- pen + pmin(1/time-etaD,eps)*pmin(1/time-etaD,eps)
        etaD <- 1/time - pmax(1/time-etaD, eps);
        logh <- etas + log(etaDs) + log(1/time -etaD)
        H <- exp(etas)
        pen - (logh*event - H)
    }
    gradi <- function(betafull) {
        beta <- betafull[1:ncol(args$X)]
        betas <- betafull[-(1:ncol(args$X))]
        time <- args$time
        eta <- as.vector(args$X %*% beta)
        etaD <- as.vector(args$XD %*% beta)
        logtstar <- log(args$time) - eta
        Xs <- predict(args$design,logtstar)
        XDs <- predict(args$designD,logtstar)
        XDDs <- predict(args$designDD,logtstar)
        etas <- as.vector(Xs %*% betas)
        etaDs <- as.vector(XDs %*% betas)
        etaDDs <- as.vector(XDDs %*% betas)
        ## H calculations
        H <- exp(etas)
        dHdbetas <- H*Xs
        dHdbeta <- -H*etaDs*X
        ## penalties
        eps <- etaDs*0. + 1e-8;
        pindexs <- etaDs < eps
        pindex <- 1/time-etaD < eps
        ## fix bounds on etaDs
        pgrads <- cbind(-2*etaDs*etaDDs*X,2*etaDs*XDs)
        etaDs <- pmax(etaDs, eps)
        ## fix bounds on etaD
        pgrad <- cbind(-2*(1/time-etaD)*XD,XDs*0)
        etaD <- 1/time - pmax(1/time-etaD, eps)
        ## 
        logh <- etas + log(etaDs) + log(1/time -etaD)
        h <- exp(logh)
        dloghdbetas <- Xs+XDs/etaDs*(!pindexs)
        dloghdbeta <- -etaDs*X*(!pindexs & !pindex) - etaDDs*X/etaDs*(!pindexs & !pindex) - XD/(1/time-etaD)*(!pindex & !pindexs)
        dhdbetas <- h*dloghdbetas
        dhdbeta <- h*dloghdbeta
        cbind(-dloghdbeta*event+dHdbeta, -dloghdbetas*event+dHdbetas) + pindex*pgrad + pindexs*pgrads
    }
    gradient2 <- function(betafull)
        colSums(gradi(betafull))
    ## browser()
    if (FALSE) {

        library(rstpm2)
        ##debug(aft)
        brcancer2 <- transform(brcancer, entry=ifelse(hormon==0, rectime/2, 0))
        system.time(aft1 <- aft(Surv(entry,rectime,censrec==1)~hormon,data=brcancer2,df=4,use.gr=FALSE))
        system.time(aft0 <- aft(Surv(entry,rectime,censrec==1)~hormon,data=brcancer2,df=4))
        
        library(rstpm2)
        system.time(aft0 <- aft(Surv(rectime,censrec==1)~hormon,data=brcancer,df=4))
        system.time(aft1 <- aft(Surv(rectime,censrec==1)~hormon,data=brcancer,df=4,use.gr=FALSE))
        ##
        brcancer100 <- brcancer[rep(1:nrow(brcancer),each=100),]
        system.time(aft0 <- aft(Surv(rectime,censrec==1)~hormon,data=brcancer100,df=4))
        system.time(aft1 <- aft(Surv(rectime,censrec==1)~hormon,data=brcancer100,df=4,use.gr=FALSE))
        ## in browser():
        args2 <- args
        args$return_type <- "nmmin"
        args2$return_type <- "vmmin"
        system.time(fit <- .Call("aft_model_output", args, PACKAGE = "rstpm2"))
        system.time(fit2 <- .Call("aft_model_output", args2, PACKAGE = "rstpm2"))
        ##
        scale <- c(1,1,1,1,1)
        as.vector(gradient(scale*init))
        colSums(gradi(scale*init))
        tmp <- sapply(1:length(init), function(i) {eps=1e-4; delta=rep(0,length(init)); delta[i]=eps; (neglli(scale*init+delta)-neglli(scale*init-delta))/(2*eps) })
        apply(tmp - gradi(scale*init), 2, range)
        ##
        ## check designD and designDD
        head(predict(design,logtstar+1e-5)-predict(design,logtstar-1e-5))/2e-5
        head(predict(designD,logtstar))
        head(predict(designD,logtstar+1e-5)-predict(designD,logtstar-1e-5))/2e-5
        head(predict(designDD,logtstar))
        ##        
        aft1 <- aft(Surv(rectime,censrec==1)~hormon,data=brcancer,df=4,
                    init=coef(aft0))
        init <- coef(aft0)
        scale <- -1
        negll(scale*init)
        negll.slow(scale*init)
        tmp <- sapply(1:length(init), function(i) {eps=1e-5; delta=rep(0,length(init)); delta[i]=eps; (neglli(scale*init+delta)-neglli(scale*init-delta))/(2*eps) })
        head(tmp[!event,])
        head(gradi(scale*init)[!event,])
        head(tmp[event,])
        head(gradi(scale*init)[event,])
        range(tmp - gradi(scale*init))
    }
    parnames(negll) <- names(init)
    ## MLE
    args$return_type <- if (use.gr) "vmmin" else "nmmin"
    fit <- .Call("aft_model_output", args, PACKAGE="rstpm2")
    args$init <- coef <- as.vector(fit$coef)
    hessian <- fit$hessian
    names(coef) <- rownames(hessian) <- colnames(hessian) <- names(init)
    mle2 <- mle2(negll, coef, gr=gradient, vecpar=TRUE, control=control, ..., eval.only=TRUE)
    mle2@vcov <- if (!inherits(vcov <- try(solve(hessian)), "try-error")) vcov else matrix(NA,length(coef), length(coef))
    mle2@details$convergence <- fit$fail # fit$itrmcd
    out <- as(mle2, "aft")
    out@args <- args
    return(out)
  }

setMethod("predict", "aft",
          function(object,newdata=NULL,
                   type=c("surv","cumhaz","hazard","density","hr","sdiff","hdiff","loghazard","link","meansurv","meansurvdiff","odds","or","meanhaz","af","fail","accfac"),
                   grid=FALSE,seqLength=300,level=0.95,
                   se.fit=FALSE,link=NULL,exposed=incrVar(var),var=NULL,keep.attributes=TRUE,...) {
              type <- match.arg(type)
              args <- object@args
              if (type %in% c("fail")) {
                  out <- 1-predict(object,newdata=newdata,type="surv",grid,seqLength,se.fit,link,exposed,var,keep.attributes,...)
                  if (se.fit) {temp <- out$lower; out$lower <- out$upper; out$upper <- temp}
                  return(out)
              }
              if (is.null(exposed) && is.null(var) & type %in% c("hr","sdiff","hdiff","meansurvdiff","or","af","accfac"))
                  stop('Either exposed or var required for type in ("hr","sdiff","hdiff","meansurvdiff","or","af","accfac")')
              ## exposed is a function that takes newdata and returns the revised newdata
              ## var is a string for a variable that defines a unit change in exposure
              if (is.null(newdata) && type %in% c("hr","sdiff","hdiff","meansurvdiff","or","af","accfac"))
                  stop("Prediction using type in ('hr','sdiff','hdiff','meansurvdiff','or','af','accfac') requires newdata to be specified.")
              calcX <- !is.null(newdata)
              time <- NULL
              if (is.null(newdata)) {
                  ##mm <- X <- model.matrix(object) # fails (missing timevar)
                  X <- args$X
                  XD <- args$XD
                  ##y <- model.response(object@model.frame)
                  y <- args$y
                  time <- as.vector(y[,ncol(y)-1])
                  newdata <- as.data.frame(args$data)
              }
              lpfunc <- function(x,...) {
                  newdata2 <- newdata
                  newdata2[[object@args$timeVar]] <- x
                  lpmatrix.lm(object@args$lm.obj,newdata2)
              }
              ## resp <- attr(Terms, "variables")[attr(Terms, "response")] 
              ## similarly for the derivatives
              if (grid) {
                  Y <- args$y
                  event <- Y[,ncol(Y)]==1 
                  time <- args$data[[args$timeVar]]
                  eventTimes <- time[event]
                  tt <- seq(min(eventTimes),max(eventTimes),length=seqLength)[-1]
                  data.x <- data.frame(tt)
                  names(data.x) <- args$timeVar
                  newdata[[args$timeVar]] <- NULL
                  newdata <- merge(newdata,data.x)
                  calcX <- TRUE
              }
              if (calcX)  {
                  X <- lpmatrix.lm(args$lm.obj, newdata)
                  XD <- grad1(lpfunc,newdata[[object@args$timeVar]],
                              log.transform=object@args$log.time.transform)
                  XD <- matrix(XD,nrow=nrow(X))
                  time <- eval(args$timeExpr,newdata)
              }
              if (type %in% c("hr","sdiff","hdiff","meansurvdiff","or","af","accfac")) {
                  newdata2 <- exposed(newdata)
                  X2 <- lpmatrix.lm(args$lm.obj, newdata2)
                  XD2 <- grad1(lpfunc,newdata2[[object@args$timeVar]],
                              log.transform=object@args$log.time.transform)
                  XD2 <- matrix(XD2,nrow=nrow(X))
                  time2 <- eval(args$timeExpr,newdata2) # is this always equal to time?
              }
              ## colMeans <- function(x) colSums(x)/apply(x,2,length)
              local <-  function (object, newdata=NULL, type="surv", exposed)
              {
                  args <- object@args
                  betafull <- coef(object)
                  beta <- betafull[1:ncol(args$X)]
                  betas <- betafull[-(1:ncol(args$X))]
                  tt <- args$terms
                  eta <- as.vector(X %*% beta)
                  logtstar <- log(time) - eta
                  etas <- as.vector(predict(args$design, logtstar) %*% betas)
                  H <- exp(etas)
                  S <- exp(-H)
                  if (type=="cumhaz") 
                      return(H)
                  if (type=="surv")
                      return(S)
                  if (type=="fail")
                      return(1-S)
                  if (type=="odds") 
                      return((1-S)/S)
                  if (type=="meansurv") 
                      return(tapply(S,newdata[[object@args$timeVar]],mean))
                  etaDs <- as.vector(predict(args$designD, logtstar) %*% betas)
                  etaD <- as.vector(XD %*% beta)
                  h <- H*etaDs*(1/time-etaD)
                  Sigma = vcov(object)
                  if (type=="link") 
                      return(eta)
                  if (type=="density")
                      return (S*h)
                  if (type=="hazard") 
                      return(h)
                  if (type=="loghazard") 
                      return(log(h))
                  if (type=="meanhaz") 
                      return(tapply(S*h,newdata[[object@args$timeVar]],sum)/tapply(S,newdata[[object@args$timeVar]],sum))
                  eta2 <- as.vector(X2 %*% beta)
                  logtstar2 <- log(time2) - eta2
                  etas2 <- as.vector(predict(args$design, logtstar2) %*% betas)
                  H2 <- exp(etas2)
                  S2 <- exp(-H2)
                  if (type=="sdiff")
                      return(S2-S)
                  if (type=="or") 
                      return((1-S2)/S2/((1-S)/S))
                  if (type=="meansurvdiff")
                      return(tapply(S2,newdata[[object@args$timeVar]],mean) - tapply(S,newdata[[object@args$timeVar]],mean))
                  etaDs2 <- as.vector(predict(args$designD, logtstar2) %*% betas)
                  etaD2 <- as.vector(XD2 %*% beta)
                  h2 <- H2*etaDs2*(1/time2-etaD2)
                  if (type=="hdiff") 
                      return(h2 - h)
                  if (type=="hr") 
                      return(h2/h)
                  if (type=="af") {
                      meanS <- tapply(S,newdata[[object@args$timeVar]],mean)
                      meanS2 <- tapply(S2,newdata[[object@args$timeVar]],mean)
                      return((meanS2 - meanS)/(1-meanS))
                  }
                  if (type=="accfac") {
                      accfac <- eta - log(1-time*etaD)
                      accfac2 <- eta2 - log(1-time2*etaD2)
                      return(exp(accfac2-accfac))
                  }
              }
              out <- if (!se.fit) {
                          local(object,newdata,type=type,exposed=exposed,
                                ...)
                      }
                      else {
                          if (is.null(link))
                              link <- switch(type,surv="cloglog",cumhaz="log",hazard="log",hr="log",sdiff="I",
                                             hdiff="I",loghazard="I",link="I",odds="log",or="log",meansurv="I",meanhaz="I",af="I",accfac="log")
                          invlinkf <- switch(link,I=I,log=exp,cloglog=cexpexp,logit=expit)
                          pred <- predictnl(object,local,link=link,newdata=newdata,type=type,gd=NULL,
                                            exposed=exposed,...)
                          ci <- confint.predictnl(pred, level = level)
                          out <- data.frame(Estimate=pred$fit,
                                            lower=ci[,1],
                                            upper=ci[,2])
                          if (link=="cloglog") 
                              out <- data.frame(Estimate=out$Estimate,lower=out$upper,upper=out$lower)
                          invlinkf(out)
                      }
              if (keep.attributes)
                  attr(out,"newdata") <- newdata
              return(out)
          })
plot.aft.meansurv <- function(x, y=NULL, times=NULL, newdata=NULL, type="meansurv", exposed=NULL, add=FALSE, ci=!add, rug=!add, recent=FALSE,
                          xlab=NULL, ylab=NULL, lty=1, line.col=1, ci.col="grey", seqLength=301, ...) {
    ## if (is.null(times)) stop("plot.meansurv: times argument should be specified")
    args <- x@args
    if (is.null(newdata)) newdata <- as.data.frame(args$data)
    if (is.null(times)) {
      Y <- args$y
      event <- Y[,ncol(Y)]==1
      time <- args$data[[args$timeVar]]
      eventTimes <- time[event]
      times <- seq(min(eventTimes),max(eventTimes),length=seqLength)[-1]
    }
    times <- times[times !=0]
    if (recent) {
        newdata <- do.call("rbind",
                           lapply(times, 
                                  function(time) {
                                      newd <- newdata
                                      newd[[args$timeVar]] <- newdata[[args$timeVar]]*0+time
                                      newd
                                  }))
        pred <- predict(x, newdata=newdata, type=type, se.fit=ci, exposed=exposed) # requires recent version
        if (type=="meansurv")
            pred <- if (ci) rbind(c(Estimate=1,lower=1,upper=1),pred) else c(1,pred)
    } else {
        pred <- lapply(times, 
                       function(time) {
                           newdata[[args$timeVar]] <- newdata[[args$timeVar]]*0+time
                           predict(x, newdata=newdata, type=type, se.fit=ci, grid=FALSE, exposed=exposed)
                       })
        pred <- do.call("rbind", pred)
        if (type=="meansurv")  {
            pred <- if (ci) rbind(c(Estimate=1,lower=1,upper=1),pred) else c(1,unlist(pred))
            times <- c(0,times)
            }
        }
    if (is.null(xlab)) xlab <- deparse(args$timeExpr)
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
        Y <- args$y
        eventTimes <- Y[Y[,ncol(Y)]==1,ncol(Y)-1]
        rug(eventTimes,col=line.col)
    }
    return(invisible(y))
}
plot.aft.base <- 
          function(x,y,newdata=NULL,type="surv",
                   xlab=NULL,ylab=NULL,line.col=1,ci.col="grey",lty=par("lty"),
                   add=FALSE,ci=!add,rug=!add,
                   var=NULL,exposed=incrVar(var),times=NULL,...) {
              if (type %in% c("meansurv","meansurvdiff","af")) {
                  return(plot.aft.meansurv(x,times=times,newdata=newdata,type=type,xlab=xlab,ylab=ylab,line.col=line.col,ci.col=ci.col,
                                       lty=lty,add=add,ci=ci,rug=rug, exposed=exposed, ...))
              }
              args <- x@args
              if (is.null(newdata)) stop("newdata argument needs to be specified")
              y <- predict(x,newdata,type=switch(type,fail="surv",margfail="margsurv",type),var=var,exposed=exposed,
                           grid=!(args$timeVar %in% names(newdata)), se.fit=ci)
              if (type %in% c("fail","margfail")) {
                  if (ci) {
                      y$Estimate <- 1-y$Estimate
                      lower <- y$lower
                      y$lower=1-y$upper
                      y$upper=1-lower
                      } else y <- structure(1-y,newdata=attr(y,"newdata"))
              }
              if (is.null(xlab)) xlab <- deparse(args$timeExpr)
              if (is.null(ylab))
                  ylab <- switch(type,hr="Hazard ratio",hazard="Hazard",surv="Survival",density="Density",
                                 sdiff="Survival difference",hdiff="Hazard difference",cumhaz="Cumulative hazard",
                                 loghazard="log(hazard)",link="Linear predictor",meansurv="Mean survival",
                                 meansurvdiff="Difference in mean survival",odds="Odds",or="Odds ratio",
                                 margsurv="Marginal survival",marghaz="Marginal hazard",marghr="Marginal hazard ratio", haz="Hazard",fail="Failure",
                                 meanhaz="Mean hazard",margfail="Marginal failure",af="Attributable fraction",meanmargsurv="Mean marginal survival",
                                 uncured="Uncured distribution","Acceleration factor")
              xx <- attr(y,"newdata")
              xx <- eval(args$timeExpr,xx) # xx[,ncol(xx)]
              if (!add) matplot(xx, y, type="n", xlab=xlab, ylab=ylab, ...)
              if (ci) {
                  polygon(c(xx,rev(xx)), c(y[,2],rev(y[,3])), col=ci.col, border=ci.col)
                  lines(xx,y[,1],col=line.col,lty=lty,...)
              } else lines(xx,y,col=line.col,lty=lty,...)
              if (rug) {
                  Y <- args$y
                  eventTimes <- Y[Y[,ncol(Y)]==1,ncol(Y)-1]
                  rug(eventTimes,col=line.col)
              }
              return(invisible(y))
          }
setMethod("plot", signature(x="aft", y="missing"),
          function(x,y,newdata=NULL,type="surv",
                   xlab=NULL,ylab=NULL,line.col=1,ci.col="grey",lty=par("lty"),
                   add=FALSE,ci=!add,rug=!add,
                   var=NULL,exposed=incrVar(var),times=NULL,...)
              plot.aft.base(x=x, y=y, newdata=newdata, type=type, xlab=xlab,
                              ylab=ylab, line.col=line.col, lty=lty, add=add,
                              ci=ci, rug=rug, var=var, exposed=exposed, times=times, ...)
          )
predictSurvival.aft <- function(obj, time=obj@args$time, X=obj@args$X) {
    localargs <- obj@args
    localargs$return_type <- "survival"
    localargs$X <- X
    localargs$time <- time
    as.vector(.Call("aft_model_output", localargs, PACKAGE="rstpm2"))
    }

## simulate from Weibull with one binary covariate
if (FALSE) {

    require(rstpm2)
    summary(aft0 <- aft(Surv(rectime,censrec==1)~hormon,data=brcancer,df=4))
    aft1 <- aft(Surv(rectime,censrec==1)~hormon,data=brcancer,df=4,init=coef(aft1))
    
    require(rstpm2)
    summary(aft0 <- aft(Surv(rectime,censrec==1)~hormon,data=brcancer,df=4))
    plot(survfit(Surv(rectime,censrec==1)~hormon,data=brcancer),col=1:2)
    plot(aft0,newdata=data.frame(hormon=0), add=TRUE, line.col="green", ci=FALSE)
    plot(aft0,newdata=data.frame(hormon=1), add=TRUE, line.col="blue", ci=FALSE)

    summary(aft1 <- aft(Surv(rectime,censrec==1)~hormon,data=brcancer,df=4,smooth.formula=~hormon:ns(log(rectime),df=3)))
    plot(survfit(Surv(rectime,censrec==1)~hormon,data=brcancer),col=1:2)
    plot(aft1,newdata=data.frame(hormon=0), add=TRUE, line.col="green", ci=FALSE)
    plot(aft1,newdata=data.frame(hormon=1), add=TRUE, line.col="blue", ci=FALSE)

    require(rstpm2)
    require(survival)
    require(splines)
    set.seed(12345)
    n <- 1e4
    x <- rep(c(0,1),n/2)
    af <- exp(0.3*x)
    time <- rweibull(n,2,5)*af
    censor <- rexp(n,rate=1/10)
    obstime <- pmin(time,censor)
    event <- ifelse(obstime==time,1,0)
    X <- cbind(x)
    XD <- X*0
    dat1 <- data.frame(obstime,event,x)
    aft1 <- aft(Surv(obstime,event)~x,data=dat1, df=2, reltol=1e-12, control=list(maxit=2000))

    plot(survfit(Surv(obstime,event)~x,data=dat1),col=1:2)
    plot(aft1,newdata=data.frame(x=0), add=TRUE, line.col="green", ci=FALSE)
    plot(aft1,newdata=data.frame(x=1), add=TRUE, line.col="blue", ci=FALSE)
    
    head(rstpm2:::predictSurvival.aft(aft1)) - head(predict(aft1))

    plot(aft1,newdata=data.frame(x=0), type="hazard", line.col="green", rug=FALSE)
    plot(aft1,newdata=data.frame(x=1), type="hazard", add=TRUE, line.col="blue", ci=TRUE)

    aft2 <- aft(Surv(obstime,event)~x,data=dat1, df=4, reltol=1e-12, control=list(maxit=2000), smooth.formula=~x:ns(log(obstime),df=3))

    plot(survfit(Surv(obstime,event)~x,data=dat1),col=1:2)
    plot(aft2,newdata=data.frame(x=0),line.col="green",add=TRUE,ci=FALSE)
    plot(aft2,newdata=data.frame(x=1),line.col="blue",add=TRUE,ci=FALSE)

    plot(aft2,newdata=data.frame(x=0),type="accfac",exposed=function(data) transform(data,x=1),se.fit=TRUE)

    aft3 <- aft(Surv(obstime,event)~x,data=dat1, df=4, reltol=1e-12, control=list(maxit=2000), smooth.formula=~x:ns(obstime,df=3))
    plot(aft3,newdata=data.frame(x=0),type="accfac",exposed=function(data) transform(data,x=1))
    plot(aft2,newdata=data.frame(x=0),type="accfac",exposed=function(data) transform(data,x=1),ci=FALSE,add=TRUE,line.col="blue")
    

    
## f$fn is not the same (?!)
    aft1$negll(aft1$par)
    -sum(aft1$logh*event-aft1$H) # ok
    -sum(aft1$f$report(aft1$par)$logh*event-aft1$f$report(aft1$par)$H) # ok
    -sum(aft2$f$report(aft1$par)$logh*event-aft2$f$report(aft1$par)$H) # ok
    aft1$f$fn(aft1$par) # ???
    
    ## f$fn is not the same (?!)
    aft2$negll(aft2$par)
    -sum(aft2$logh*event-aft2$H) # ok
    -sum(aft2$f$report(aft2$par)$logh*dat1$event-aft2$f$report(aft2$par)$H) # ok
    aft2$f$fn(aft2$par) # ???

    ## the events are the same
    all(event == aft1$f$report(aft1$par)$event) # ok
    all(event == aft2$f$report(aft2$par)$event) # ok
    ## The H and logh vectors are very close
    max(abs(aft1$f$report(aft1$par)$H - aft1$H))       # ok
    max(abs(aft2$f$report(aft2$par)$H - aft2$H))       # ok
    max(abs(aft1$f$report(aft1$par)$logh - aft1$logh)) # ok
    max(abs(aft2$f$report(aft2$par)$logh - aft2$logh)) # ok
    ##
    ## the Xs and XDs matrices are very close
    max(abs(aft1$f$report(aft1$par)$Xs - aft1$Xs))   # ok
    max(abs(aft1$f$report(aft1$par)$XDs - aft1$XDs)) # ok
    max(abs(aft2$f$report(aft2$par)$Xs - aft2$Xs))   # ok
    max(abs(aft2$f$report(aft2$par)$XDs - aft2$XDs)) # ok


    head(aft1$Xs)
    head(aft2$Xs)
    head(aft1$f$report(aft1$par)$Xs)
    head(aft1$f$report(aft2$par)$Xs)
    head(aft2$f$report(aft2$par)$Xs)
    
    aft1$f$gr(aft1$par) # ok
    delta=function(i,eps=1e-5) {new=rep(0,5); new[i]=eps; new}
    for (i in 1:length(aft1$par))
        print((aft1$f$fn(aft1$par+delta(i))-aft1$f$fn(aft1$par-delta(i)))/(2e-5))
    ## Gradient at the 'final' value are not the same
    for (i in 1:length(aft2$par))
        print((aft2$negll(aft1$par+delta(i))-aft2$negll(aft1$par-delta(i)))/(2e-5))

    ## gradient at the initial value are the same
    aft1$f$gr(aft1$init) # ok
    delta=function(i,eps=1e-5) {new=rep(0,5); new[i]=eps; new}
    for (i in 1:length(aft1$par))
        print((aft1$f$fn(aft1$init+delta(i))-aft1$f$fn(aft1$init-delta(i)))/(2e-5))

    ## Objective at the initial values are the same
    as.numeric(aft1$negll(aft1$init)) - aft1$f$fn(aft1$init) # ok
    as.numeric(aft2$negll(aft2$init)) - aft2$f$fn(aft2$init) # ok

    ## objectives at the 'final' values are NOT the same
    as.numeric(aft1$negll(aft1$init+0.1)) - aft1$f$fn(aft1$init+0.1) # ??
    as.numeric(aft1$negll(aft1$par)) - aft1$f$fn(aft1$par) # ???

    undebug(aft1$negll)
    aft1$negll(aft1$par)

    aft1$init
    aft1$f$par
    aft1$f$fn(aft1$init)
    aft1$f$fn(aft1$f$par)
    aft1$f$fn(aft1$par)
    aft1$f$fn(aft2$par)
    
}
