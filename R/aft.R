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
                tvc = NULL, cure.formula=~1,
                control = list(), init = NULL,
                weights = NULL, tvc.intercept=TRUE, tvc.integrated= FALSE,
                timeVar = "", time0Var = "",
                cure = FALSE, mixture = FALSE,
                contrasts = NULL, subset = NULL, ...) {
    dots = list(...)
    control.defaults = list(parscale = 1, maxit = 1000, constrOptim = FALSE,
                            use.gr = TRUE, reltol = 1.0e-8, trace = 0,
                            nNodes = 20, add.penalties = TRUE)
    control = modifyList(control.defaults, control)
    if (any(ioverlap <- names(dots) %in% names(control.defaults))) {
        overlap = names(dots)[ioverlap]
        for (name in overlap) {
            cat(sprintf("Deprecated argument: %s; use control argument\n", name))
            control[[name]] = dots[[name]]
        }
    }
    ## Special case
    if (mixture && df>2) {
        Call = match.call()
        Call$df = 2
        fitWeibullMixture = eval(Call,parent.frame())
    } else {
        fitWeibullMixture = NULL
    }
    ## parse the event expression
    eventInstance <- eval(lhs(formula),envir=data)
    stopifnot(length(lhs(formula))>=2)
    eventExpr <- lhs(formula)[[length(lhs(formula))]]
    delayed <- length(lhs(formula))>=4 # indicator for multiple times (cf. strictly delayed)
    surv.type <- attr(eventInstance,"type")
    if (surv.type %in% c("interval","interval2","left","mstate"))
        stop("aft not implemented for Surv type ",surv.type,".")
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
    if (!is.null(tvc)) {
        tvc.formulas <-
            lapply(names(tvc), function(name)
                call(":",
                     call("as.numeric",as.name(name)),
                     as.call(c(quote(ns),
                               if (tvc.integrated) timeExpr else call("log",timeExpr),
                               vector2call(list(intercept=tvc.intercept,df=tvc[[name]]))))))
        if (length(tvc.formulas)>1)
            tvc.formulas <- list(Reduce(`%call+%`, tvc.formulas))
        tvc.formula <- as.formula(call("~",tvc.formulas[[1]]))
        rhs(full.formula) <- rhs(full.formula) %call+% rhs(tvc.formula)
    }
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
    ## Specials:
    bhazard = rep(0,nrow(data))
    specials.names <- c("bhazard")
    specials <- attr(terms.formula(full.formula, specials.names), "specials")
    spcall <- mf
    spcall[[1]] <- quote(stats::model.frame)
    spcall$formula <- terms(formula, specials.names, data = data)
    mf2 <- eval(spcall, parent.frame())
    lm.formula = full.formula
    if (any(!sapply(specials,is.null))) {
        bhazard.index <- specials$bhazard
        if (length(bhazard.index)>0) {
            bhazard <- mf2[, bhazard.index]
            bhazard.index2 <- attr(terms.formula(full.formula, "bhazard"), "specials")$bhazard
            termobj = terms(mf2)
            dropped.terms = if(length(attr(termobj,"term.labels"))==1) reformulate("1", response=termobj[[2L]], intercept=TRUE, env=environment(termobj)) else stats::drop.terms(terms(mf2), bhazard.index - 1, keep.response=TRUE)
            formula <- formula(dropped.terms)
            lm.formula <- formula(stats::drop.terms(terms(full.formula), bhazard.index2 - 1))
        } else {
            bhazard.index2 = NULL
        }
        ## rm(mf2,spcall)
    }
    ## get variables
    time <- eval(timeExpr, data, parent.frame())
    time0Expr <- NULL # initialise
    if (delayed) {
        time0Expr <- lhs(formula)[[2]]
        if (time0Var == "")
            time0Var <- all.vars(time0Expr)
        time0 <- eval(time0Expr, data, parent.frame())
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
    coxph.call$formula = formula
    coxph.call$model <- TRUE
    coxph.obj <- eval(coxph.call, envir=parent.frame())
    y <- model.extract(model.frame(coxph.obj),"response")
    data$logHhat <- pmax(-18,log(-log(S0hat(coxph.obj))))
    ## now for the cure fraction
    glm.cure.call = coxph.call
    glm.cure.call[[1]] = as.name("glm")
    glm.cure.call$family = as.name("binomial")
    glm.cure.call$data = as.name("data")
    glm.cure.call$model = NULL
    data["*event*"] = event 
    lhs(glm.cure.call$formula) = as.name("*event*")
    rhs(glm.cure.call$formula) = rhs(cure.formula)
    ## glm(y ~ X, family=binomial)
    glm.cure.obj <- eval(glm.cure.call, envir=as.list(environment()), enclos=parent.frame())
    Xc = model.matrix(glm.cure.obj, data)
    ##
    ## pred1 <- predict(survreg1)
    data$logtstar <- log(time)
    ## data$logtstar <- log(time/pred1)
    ## initial values and object for lpmatrix predictions
    lm.call <- mf
    lm.call[[1L]] <- as.name("lm")
    ## lm.formula <- full.formula # ??
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
    ## lm0D.obj <- lm(logHhat~nsxD(logtstar,df,intercept=TRUE,cure=cure)-1,dataEvents)
    coef0 <- coef(lm0.obj) # log-log baseline
    ## design information for baseline survival
    design <- nsx(dataEvents$logtstar, df=df, intercept=TRUE, cure=cure)
    designD <- nsxD(dataEvents$logtstar, df=df, intercept=TRUE, cure=cure)
    designDD <- nsxDD(dataEvents$logtstar, df=df, intercept=TRUE, cure=cure)
    ##
    ## set up mf and wt
    mt <- terms(lm.obj)
    mf <- model.frame(lm.obj)
    ## wt <- model.weights(lm.obj$model)
    wt <- if (is.null(substitute(weights))) rep(1,nrow(data)) else eval(substitute(weights),data,parent.frame())
    ##
    ## XD matrix
    loglpfunc <- function(x,fit,data,var) {
        data[[var]] <- exp(x)
        lpmatrix.lm(fit,data)
    }
    ##
    X <- lpmatrix.lm(lm.obj,data)
    if (is.integer(X.index <- which.dim(X)))
        warning("Design matrix for the acceleration factor is not full rank")
    X <- X[, X.index, drop=FALSE]
    XD0 <- X0 <- XD <- matrix(0,1,ncol(X))
    X_list = list()
    gauss = gauss.quad(control$nNodes)
    if (delayed && all(time0==0)) delayed <- FALSE # CAREFUL HERE: delayed redefined
    if (tvc.integrated) {
        X_list = lapply(1:control$nNodes, function(i)
            lpmatrix.lm(lm.obj,
                        local({ data[[timeVar]] = (gauss$nodes[i]+1)/2*data[[timeVar]]; data}))[,X.index, drop = FALSE])
        if (delayed) {
            X_list0 = lapply(1:control$nNodes, function(i)
                lpmatrix.lm(lm.obj,
                            local({ data[[timeVar]] = (gauss$nodes[i]+1)/2*data[[time0Var]]; data}))[,X.index, drop = FALSE])
        }
    } else { # tvc using cumulative formulation 
        XD <- grad1(loglpfunc,log(data[[timeVar]]),lm.obj,data,timeVar,log.transform=FALSE)
        XD <- matrix(XD,nrow=nrow(X))[,X.index, drop = FALSE]
        if (delayed) {
            ind0 <- time0>0
            data0 <- data[ind0,,drop=FALSE] # data for delayed entry times
            data0[[timeVar]] <- data0[[time0Var]]
            X0 <- lpmatrix.lm(lm.obj, data0)
            XD0 <- grad1(loglpfunc,log(data0[[timeVar]]),lm.obj,data0,timeVar,
                         log.transform=FALSE)
            XD0 <- matrix(XD0,nrow=nrow(X))[,X.index, drop = FALSE]
            rm(data0)
        }
    }
    ## Weibull regression
    if (delayed) {
        if (requireNamespace("eha", quietly = TRUE)) {
            survreg1 <- eha::aftreg(formula, data)
            coef1 <- -coef(survreg1) # reversed parameterisation
            coef1 <- coef1[1:(length(coef1)-2)][X.index]
        } else coef1 <- rep(0,ncol(X))
    } else {
        survreg1 <- survival::survreg(formula, data)
        coef1 <- coef(survreg1)
        coef1 <- coef1[-1] # assumes intercept included in the formula
    }
    if (ncol(X)>length(coef1)) {
        coef1 <- c(coef1,rep(0,ncol(X) - length(coef1)))
        names(coef1) <- names(coef1b)[X.index]
    }
    coef2 = coef(glm.cure.obj)
    names(coef2) = paste0("cure.", names(coef2))
    if (is.null(init))
        init <- if (mixture) c(coef1, -coef2, coef0) else c(coef1, coef0) # -coef2 because the glm models for uncured!
    if (any(is.na(init) | is.nan(init)))
        stop("Some missing initial values - check that the design matrix is full rank.")
    if (!is.null(control) && "parscale" %in% names(control)) {
        if (length(control$parscale)==1)
            control$parscale <- rep(control$parscale,length(init))
        if (is.null(names(control$parscale)))
            names(control$parscale) <- names(init)
    }
    parscale <- rep(if (is.null(control$parscale)) 1 else control$parscale,length=length(init))
    names(parscale) <- names(init)
    args <- list(init=init,X=X,XD=XD,
                 X_list=X_list,
                 X_list0=if (delayed) X_list0 else list(matrix(0,0,0)),
                 X.index=X.index,
                 wt=wt,event=ifelse(event,1,0),time=time,y=y,
                 time0 = if (delayed) time0 else 0*time,
                 timeVar=timeVar,timeExpr=timeExpr,terms=mt,
                 delayed=delayed, X0=X0, XD0=XD0,
                 parscale=parscale, reltol=control$reltol,
                 Xc=if (mixture) Xc else matrix(0,0,0), maxit=control$maxit,
                 trace = as.integer(control$trace),
                 boundaryKnots=attr(design,"Boundary.knots"), q.const=t(attr(design,"q.const")),
                 interiorKnots=attr(design,"knots"), design=design, designD=designD,
                 designDD=designDD, cure=as.integer(cure), mixture = as.integer(mixture),
                 tvc.integrated=tvc.integrated,
                 data=data, lm.obj = lm.obj, glm.cure.obj = glm.cure.obj,
                 init_copy = init, return_type="optim",
                 gweights=gauss$weights, gnodes=gauss$nodes, bhazard=bhazard,
                 add.penalties = control$add.penalties, df=df, ci=rep(0, df+1),
                 constrOptim = control$constrOptim)
    getui = function(args) {
        q.const = t(args$q.const) # (df+2) x df
        n.coef = length(args$init)
        nr = nrow(q.const)
        df = ncol(q.const)
        (cbind(0,diag(nr-1))-cbind(diag(nr-1),0)) %*%
            cbind(matrix(0,nr,n.coef-df), q.const) # (df+1) x n.coef
    }
    args$ui = getui(args)
    negll <- function(beta, ...) {
        stopifnot(is.numeric(beta),
                  length(beta) == length(args$init))
        localargs <- args
        localargs$return_type <- "objective"
        localargs$init <- beta
        if (length(dots <- list(...)) > 0)
            localargs <- modifyList(localargs, dots)
        return(.Call("aft_model_output", localargs, PACKAGE="rstpm2"))
    }
    gradient <- function(beta, ...) {
        stopifnot(is.numeric(beta),
                  length(beta) == length(args$init))
        localargs <- args
        localargs$return_type <- "gradient"
        localargs$init <- beta
        if (length(dots <- list(...)) > 0)
            localargs <- modifyList(localargs, dots)
        return(as.vector(.Call("aft_model_output", localargs, PACKAGE="rstpm2")))
    }
    parnames(negll) <- names(init)
    args$negll = negll
    args$gradient = gradient
    if (!is.null(fitWeibullMixture)) {
        this.index = (length(coef1)+1):(length(coef1)+length(coef2))
        fixed = as.list(coef(fitWeibullMixture)[this.index])
        this.control = control
        for (name in c("constrOptim", "use.gr", "nNodes", "add.penalties"))
            this.control[[name]] = NULL
        this.control$parscale = this.control$parscale[-this.index]
        fixedfit <- if (control$use.gr) {
                    bbmle::mle2(negll, init, vecpar=TRUE, control=this.control,
                                fixed=fixed, gr=gradient, ...)
                } else {
                    bbmle::mle2(negll, init, vecpar=TRUE, control=this.control,
                                fixed=fixed, ...)
                }
        init = coef(fixedfit)
        rm(fixedfit)
    }
    ## MLE
    if (delayed && control$use.gr) { # initial search using nmmin (conservative -- is this needed?)
        args$return_type <- "nmmin"
        args$maxit <- 50
        fit <- .Call("aft_model_output", args, PACKAGE="rstpm2")
        args$maxit <- control$maxit
    }
    optim_step <- function(use.gr) {
        args$return_type <<- if (use.gr) "vmmin" else "nmmin"
        fit <- .Call("aft_model_output", args, PACKAGE="rstpm2")
        coef <- as.vector(fit$coef)
        hessian <- fit$hessian
        names(coef) <- rownames(hessian) <- colnames(hessian) <- names(init)
        args$init <<- coef
        ## we could use mle2() to calculate vcov by setting eval.only=FALSE
        mle2 <- if (use.gr) bbmle::mle2(negll, coef, vecpar=TRUE, control=control,
                                        gr=gradient, ..., eval.only=TRUE)
                else bbmle::mle2(negll, coef, vecpar=TRUE, control=control, ..., eval.only=TRUE)
        mle2@details$convergence <- fit$fail # fit$itrmcd
        vcov <- try(solve(hessian,tol=0), silent=TRUE)
        if (inherits(vcov, "try-error"))
            vcov <- try(solve(hessian+1e-6*diag(nrow(hessian)), tol=0), silent=TRUE)
        if (inherits(vcov, "try-error")) {
            if (!use.gr)
                message("Non-invertible Hessian")
            mle2@vcov <- matrix(NA,length(coef), length(coef))
        } else {
            mle2@vcov <- vcov
        }
        mle2
    }
    ## mle2 <-  bbmle::mle2(negll, init, vecpar=TRUE, control=control, ...)
    mle2 <- optim_step(control$use.gr)
    if (all(is.na(mle2@vcov)) && control$use.gr) {
        args$init <- init
        mle2 <- optim_step(FALSE)
    }
    out <- as(mle2, "aft")
    out@args <- args
    attr(out,"nobs") <- length(out@args$event) # for logLik method
    return(out)
}

setMethod("nobs", "aft", function(object, ...) length(object@args$event))

setMethod("predict", "aft",
          function(object,newdata=NULL,
                   type=c("surv","cumhaz","hazard","density","hr","sdiff","hdiff","loghazard","link","meansurv","meansurvdiff","odds","or","meanhaz","af","fail","accfac","gradh"),
                   grid=FALSE,seqLength=300,level=0.95,
                   se.fit=FALSE,link=NULL,exposed=incrVar(var),var=NULL,keep.attributes=TRUE,...) {
              type = match.arg(type)
              if (object@args$tvc.integrated)
                  predict.aft_integrated2(object, newdata, type, grid, seqLength, level, se.fit, link, exposed, var, keep.attributes, ...)
              else predict.aft_mixture2(object, newdata, type, grid, seqLength, level, se.fit, link, exposed, var, keep.attributes, ...)
          })

predict.aft_mixture2 <-
    function(object,newdata=NULL,
             type=c("surv","cumhaz","hazard","density","hr","sdiff","hdiff","loghazard","link",
                    "meansurv","meansurvdiff","odds","or","meanhaz","af","fail","accfac","gradh"),
             grid=FALSE,seqLength=300,level=0.95,
             se.fit=FALSE,link=NULL,exposed=incrVar(var),var=NULL,keep.attributes=TRUE,...) {
        type <- match.arg(type)
        args <- object@args
        if (type %in% c("fail")) {
            out <- 1-predict(object,newdata=newdata,type="surv",grid,seqLength,se.fit,link,
                             exposed,var,keep.attributes,...)
            if (se.fit) {temp <- out$lower; out$lower <- out$upper; out$upper <- temp}
            return(out)
        }
        if (is.null(exposed) && is.null(var) & type %in% c("hr","sdiff","hdiff","meansurvdiff",
                                                           "or","af","accfac"))
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
        loglpfunc <- function(x,newdata,...) {
            newdata[[object@args$timeVar]] <- exp(x)
            lpmatrix.lm(object@args$lm.obj,newdata=newdata)
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
            X <- lpmatrix.lm(args$lm.obj, newdata)[,args$X.index, drop = FALSE]
            XD <- grad1(loglpfunc,log(newdata[[object@args$timeVar]]),
                        log.transform=FALSE, newdata=newdata)
            XD <- matrix(XD,nrow=nrow(X))[,args$X.index, drop = FALSE]
            Xc <- lpmatrix.lm(args$glm.cure.obj, newdata)
            time <- eval(args$timeExpr,newdata)
        }
        if (type %in% c("hr","sdiff","hdiff","meansurvdiff","or","af","accfac")) {
            newdata2 <- exposed(newdata)
            X2 <- lpmatrix.lm(args$lm.obj, newdata2)[,args$X.index, drop = FALSE]
            XD2 <- grad1(loglpfunc,log(newdata2[[object@args$timeVar]]),
                         log.transform=FALSE, newdata=newdata2)
            XD2 <- matrix(XD2,nrow=nrow(X))[,args$X.index, drop = FALSE]
            time2 <- eval(args$timeExpr,newdata2) # is this always equal to time?
            Xc2 = model.matrix(args$glm.cure.obj, newdata2)
        }
        if (type %in% c("grad")) {
            return(predict.aft.ext(object, type=type, time=time, X=X, XD=XD))
        }
        ## colMeans <- function(x) colSums(x)/apply(x,2,length)
        local <-  function (object, newdata=NULL, type="surv", exposed, ...)
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
                accfac <- -eta + log(1-etaD)
                accfac2 <- -eta2 + log(1-etaD2)
                return(exp(accfac2-accfac))
            }
        }
        local2 <-  function (object, newdata=NULL, type="surv", exposed, ...)
        {
            args <- object@args
            betafull <- coef(object)
            beta <- betafull[1:ncol(args$X)]
            betac <- betafull[(ncol(args$X)+1):(ncol(args$X)+ncol(args$Xc))]
            betas <- betafull[-(1:(ncol(args$X)+ncol(args$Xc)))]
            tt <- args$terms
            eta <- as.vector(X %*% beta)
            etac <- as.vector(Xc %*% betac)
            cure_frac <- exp(etac)/(1+exp(etac))
            logtstar <- log(time) - eta
            etas <- as.vector(predict(args$design, logtstar) %*% betas)
            Hu <- exp(etas)
            Su <- exp(-Hu)
            S <- cure_frac + (1-cure_frac)*Su
            if (type=="cumhaz")
                return(-log(cure_frac + (1-cure_frac)*exp(-Hu)))
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
            hu <- Hu*etaDs*(1/time-etaD)
            h <- (1-cure_frac)*exp(-Hu)*hu/(cure_frac + (1-cure_frac)*exp(-Hu))
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
            etac2 <- as.vector(Xc2 %*% betac)
            cure_frac2 <- exp(etac2)/(1+exp(etac2))
            Hu2 <- exp(etas2)
            Su2 <- exp(-Hu2)
            S2 <- cure_frac2 + (1-cure_frac2)*Su2
            if (type=="sdiff")
                return(S2-S)
            if (type=="or")
                return((1-S2)/S2/((1-S)/S))
            if (type=="meansurvdiff")
                return(tapply(S2,newdata[[object@args$timeVar]],mean) - tapply(S,newdata[[object@args$timeVar]],mean))
            etaDs2 <- as.vector(predict(args$designD, logtstar2) %*% betas)
            etaD2 <- as.vector(XD2 %*% beta)
            hu2 <- Hu2*etaDs2*(1/time2-etaD2)
            h2 <- (1-cure_frac2)*exp(-Hu2)*hu2/(cure_frac2 + (1-cure_frac2)*exp(-Hu2))
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
                accfac <- -eta + log(1-etaD)
                accfac2 <- -eta2 + log(1-etaD2)
                return(exp(accfac2-accfac))
            }
        }
        local <- if (args$mixture) local2 else local
        out <- if (!se.fit) {
                   local(object,newdata,type=type,exposed=exposed,
                         ...)
               } else {
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
    }

predict.aft_integrated2 =  function(object,newdata=NULL,
                                type=c("surv","cumhaz","hazard","density","hr","sdiff","hdiff","loghazard","link","meansurv","meansurvdiff","odds","or","meanhaz","af","fail","accfac","gradh"),
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
        X_list <- args$X_list
        XD <- args$XD
        ##y <- model.response(object@model.frame)
        y <- args$y
        time <- as.vector(y[,ncol(y)-1])
        newdata <- as.data.frame(args$data)
    }
    loglpfunc <- function(x,newdata,...) {
        newdata[[object@args$timeVar]] <- exp(x)
        lpmatrix.lm(object@args$lm.obj,newdata=newdata)
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
        X <- lpmatrix.lm(args$lm.obj, newdata)[,args$X.index, drop = FALSE]
        XD <- grad1(loglpfunc,log(newdata[[object@args$timeVar]]),
                    log.transform=FALSE, newdata=newdata)
        XD <- matrix(XD,nrow=nrow(X))[,args$X.index, drop = FALSE]
        Xc <- lpmatrix.lm(args$glm.cure.obj, newdata)
        time <- eval(args$timeExpr,newdata)
        X_list = lapply(args$gnodes, function(gnode)
            lpmatrix.lm(args$lm.obj,
                    local({ newdata[[args$timeVar]] = (gnode+1)/2*newdata[[args$timeVar]]; newdata}))[,args$X.index, drop=FALSE])
    }
    if (type %in% c("hr","sdiff","hdiff","meansurvdiff","or","af","accfac")) {
        newdata2 <- exposed(newdata)
        X2 <- lpmatrix.lm(args$lm.obj, newdata2)[,args$X.index, drop = FALSE]
        XD2 <- grad1(loglpfunc,log(newdata2[[object@args$timeVar]]),
                     log.transform=FALSE, newdata=newdata2)
        XD2 <- matrix(XD2,nrow=nrow(X))[,args$X.index, drop = FALSE]
        time2 <- eval(args$timeExpr,newdata2) # is this always equal to time?
        Xc2 = model.matrix(args$glm.cure.obj, newdata2)
        X_list2 = lapply(args$gnodes, function(gnode)
            lpmatrix.lm(args$lm.obj,
                    local({ newdata2[[args$timeVar]] = (gnode+1)/2*newdata2[[args$timeVar]]; newdata2}))[,args$X.index, drop=FALSE])
    }
    if (type == "gradh") {
        return(predict.aft.ext(object, type="gradh", time=time, X=X, XD=XD,
                                          Xc=Xc, X_list=X_list))
    }
    ## colMeans <- function(x) colSums(x)/apply(x,2,length)
    local <-  function (object, newdata=NULL, type="surv", exposed, ...)
    {
        args <- object@args
        betafull <- coef(object)
        beta <- betafull[1:ncol(X)]
        betas <- betafull[-(1:ncol(X))]
        tt <- args$terms
        tstar = time*0
        scale = time/2.0
        for(i in 1:length(X_list)) {
            tstar = tstar + args$gweights[i]*scale * exp(-X_list[[i]] %*% beta)
        }
        ## eta <- as.vector(X %*% beta)
        logtstar <- log(tstar)
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
        h <- drop(H*etaDs*exp(X %*% beta))
        Sigma = vcov(object)
        if (type=="link")
            return(logtstar) ## what should this be?
        if (type=="density")
            return (S*h)
        if (type=="hazard")
            return(h)
        if (type=="loghazard")
            return(log(h))
        if (type=="meanhaz")
            return(tapply(S*h,newdata[[object@args$timeVar]],sum)/tapply(S,newdata[[object@args$timeVar]],sum))
        tstar2 = time2*0
        scale2 = time2/2.0
        for(i in 1:length(X_list2)) {
            tstar2 = tstar2 + args$gweights[i]*scale2 * exp(-X_list2[[i]] %*% beta)
        }
        logtstar2 = log(tstar2)
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
        h2 <- drop(H2*etaDs2*exp(X2 %*% beta))
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
            return(exp(-(X2-X) %*% beta))
        }
    }
    local2 <-  function (object, newdata=NULL, type="surv", exposed, ...)
    {
        args <- object@args
        betafull <- coef(object)
        beta <- betafull[1:ncol(args$X)]
        betac <- betafull[(ncol(args$X)+1):(ncol(args$X)+ncol(args$Xc))]
        betas <- betafull[-(1:(ncol(args$X)+ncol(args$Xc)))]
        tt <- args$terms
        tt <- args$terms
        tstar = time*0
        scale = time/2.0
        for(i in 1:length(X_list)) {
            tstar = tstar + args$gweights[i]*scale * exp(-X_list[[i]] %*% beta)
        }
        ## eta <- as.vector(X %*% beta)
        etac <- as.vector(Xc %*% betac)
        cure_frac <- exp(etac)/(1+exp(etac))
        logtstar <- log(tstar)
        etas <- as.vector(predict(args$design, logtstar) %*% betas)
        Hu <- exp(etas)
        Su <- exp(-Hu)
        S <- cure_frac + (1-cure_frac)*Su
        if (type=="cumhaz")
            return(-log(cure_frac + (1-cure_frac)*exp(-Hu)))
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
        hu <- drop(Hu*etaDs*exp(X %*% beta))
        h <- (1-cure_frac)*exp(-Hu)*hu/(cure_frac + (1-cure_frac)*exp(-Hu))
        Sigma = vcov(object)
        if (type=="link")
            return(logtstar) # is this correct?
        if (type=="density")
            return (S*h)
        if (type=="hazard")
            return(h)
        if (type=="loghazard")
            return(log(h))
        if (type=="meanhaz")
            return(tapply(S*h,newdata[[object@args$timeVar]],sum)/tapply(S,newdata[[object@args$timeVar]],sum))
        tstar2 = time2*0
        scale2 = time2/2.0
        for(i in 1:length(X_list2)) {
            tstar2 = tstar2 + args$gweights[i]*scale2 * exp(-X_list2[[i]] %*% beta)
        }
        logtstar2 = log(tstar2)
        ## eta2 <- as.vector(X2 %*% beta)
        etas2 <- as.vector(predict(args$design, logtstar2) %*% betas)
        etac2 <- as.vector(Xc2 %*% betac)
        cure_frac2 <- exp(etac2)/(1+exp(etac2))
        Hu2 <- exp(etas2)
        Su2 <- exp(-Hu2)
        S2 <- cure_frac2 + (1-cure_frac2)*Su2
        if (type=="sdiff")
            return(S2-S)
        if (type=="or")
            return((1-S2)/S2/((1-S)/S))
        if (type=="meansurvdiff")
            return(tapply(S2,newdata[[object@args$timeVar]],mean) - tapply(S,newdata[[object@args$timeVar]],mean))
        etaDs2 <- as.vector(predict(args$designD, logtstar2) %*% betas)
        etaD2 <- as.vector(XD2 %*% beta)
        hu2 <- drop(Hu2*etaDs2*exp(X2 %*% beta))
        h2 <- (1-cure_frac2)*exp(-Hu2)*hu2/(cure_frac2 + (1-cure_frac2)*exp(-Hu2))
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
            return(exp(-(X2-X) %*% beta))
        }
    }
    local <- if(args$mixture)  local2 else local
    out <- if (!se.fit) {
               local(object,newdata,type=type,exposed=exposed,
                     ...)
           } else {
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
}

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

predict.aft.ext <- function(obj, type=c("survival","haz","gradh"),
                            time=obj@args$time, X=obj@args$X, XD=obj@args$XD,
                            X_list=obj@args$X_list, Xc=obj@args$Xc) {
    type <- match.arg(type)
    localargs <- obj@args
    localargs$return_type <- type
    localargs$X <- X
    localargs$XD <- XD
    localargs$time <- time
    localargs$Xc <- Xc
    localargs$X_list <- X_list
    as.matrix(.Call("aft_model_output", localargs, PACKAGE="rstpm2"))
}

## ## simulate from Weibull with one binary covariate
## if (FALSE) {
##     require(rstpm2)
##     summary(aft0 <- aft(Surv(rectime,censrec==1)~hormon,data=brcancer,df=4))
##     aft1 <- aft(Surv(rectime,censrec==1)~hormon,data=brcancer,df=4,init=coef(aft1))
##     ##
##     require(rstpm2)
##     summary(aft0 <- aft(Surv(rectime,censrec==1)~hormon,data=brcancer,df=4))
##     plot(survfit(Surv(rectime,censrec==1)~hormon,data=brcancer),col=1:2)
##     plot(aft0,newdata=data.frame(hormon=0), add=TRUE, line.col="green", ci=FALSE)
##     plot(aft0,newdata=data.frame(hormon=1), add=TRUE, line.col="blue", ci=FALSE)
##     ##
##     summary(aft1 <- aft(Surv(rectime,censrec==1)~hormon,data=brcancer,df=4,smooth.formula=~hormon:ns(log(rectime),df=3)))
##     plot(survfit(Surv(rectime,censrec==1)~hormon,data=brcancer),col=1:2)
##     plot(aft1,newdata=data.frame(hormon=0), add=TRUE, line.col="green", ci=FALSE)
##     plot(aft1,newdata=data.frame(hormon=1), add=TRUE, line.col="blue", ci=FALSE)

##     require(rstpm2)
##     require(survival)
##     require(splines)
##     set.seed(12345)
##     n <- 1e4
##     x <- rep(c(0,1),n/2)
##     af <- exp(0.3*x)
##     time <- rweibull(n,2,5)*af
##     censor <- rexp(n,rate=1/10)
##     obstime <- pmin(time,censor)
##     event <- ifelse(obstime==time,1,0)
##     X <- cbind(x)
##     XD <- X*0
##     dat1 <- data.frame(obstime,event,x)
##     aft1 <- aft(Surv(obstime,event)~x,data=dat1, df=2, reltol=1e-12, control=list(maxit=2000))

##     plot(survfit(Surv(obstime,event)~x,data=dat1),col=1:2)
##     plot(aft1,newdata=data.frame(x=0), add=TRUE, line.col="green", ci=FALSE)
##     plot(aft1,newdata=data.frame(x=1), add=TRUE, line.col="blue", ci=FALSE)
    
##     head(rstpm2:::predict.aft.ext(aft1) - predict(aft1))
##     range(rstpm2:::predict.aft.ext(aft1,type="haz") - predict(aft1, type="haz"))
##     rstpm2:::predict.aft.ext(aft1,type="haz",time=aft1@args$time[1:6],X=aft1@args$X[1:6,,drop=FALSE],XD=aft1@args$XD[1:6,,drop=FALSE]) - head(predict(aft1, type="haz"))
##     predict.aft.ext.test <- function(obj, eps=1e-5) {
##         localargs <- obj@args
##         localargs$return_type <- "haz"
##         basecoef <- coef(obj)
##         sapply(1:length(basecoef),
##                function(i) {
##                    coef <- basecoef
##                    coef[i] <- coef[i]+eps
##                    localargs$init <- coef
##                    upper <- as.vector(.Call("aft_model_output", localargs, PACKAGE="rstpm2"))
##                    coef <- basecoef
##                    coef[i] <- coef[i]-eps
##                    localargs$init <- coef
##                    lower <- as.vector(.Call("aft_model_output", localargs, PACKAGE="rstpm2"))
##                    (upper-lower)/eps/2
##                })
##     }
##     temp <- predict.aft.ext.test(aft1)
##     range(rstpm2:::predict.aft.ext(aft1,type="gradh") - temp)
##     rstpm2:::predict.aft.ext(aft1,type="gradh") - head(temp)
##     range(predict(aft1,newdata=data.frame(obstime=aft1@args$time[1:6],x=x[1:6]),type="gradh") -
##           head(temp))
    
##     plot(aft1,newdata=data.frame(x=0), type="hazard", line.col="green", rug=FALSE)
##     plot(aft1,newdata=data.frame(x=1), type="hazard", add=TRUE, line.col="blue", ci=TRUE)

##     aft2 <- aft(Surv(obstime,event)~x,data=dat1, df=4, reltol=1e-12, control=list(maxit=2000), smooth.formula=~x:ns(log(obstime),df=3))

##     plot(survfit(Surv(obstime,event)~x,data=dat1),col=1:2)
##     plot(aft2,newdata=data.frame(x=0),line.col="green",add=TRUE,ci=FALSE)
##     plot(aft2,newdata=data.frame(x=1),line.col="blue",add=TRUE,ci=FALSE)

##     plot(aft2,newdata=data.frame(x=0),type="accfac",exposed=function(data) transform(data,x=1),se.fit=TRUE)

##     aft3 <- aft(Surv(obstime,event)~x,data=dat1, df=4, reltol=1e-12, control=list(maxit=2000), smooth.formula=~x:ns(obstime,df=3))
##     plot(aft3,newdata=data.frame(x=0),type="accfac",exposed=function(data) transform(data,x=1))
##     plot(aft2,newdata=data.frame(x=0),type="accfac",exposed=function(data) transform(data,x=1),ci=FALSE,add=TRUE,line.col="blue")
    

    
## ## f$fn is not the same (?!)
##     aft1$negll(aft1$par)
##     -sum(aft1$logh*event-aft1$H) # ok
##     -sum(aft1$f$report(aft1$par)$logh*event-aft1$f$report(aft1$par)$H) # ok
##     -sum(aft2$f$report(aft1$par)$logh*event-aft2$f$report(aft1$par)$H) # ok
##     aft1$f$fn(aft1$par) # ???
    
##     ## f$fn is not the same (?!)
##     aft2$negll(aft2$par)
##     -sum(aft2$logh*event-aft2$H) # ok
##     -sum(aft2$f$report(aft2$par)$logh*dat1$event-aft2$f$report(aft2$par)$H) # ok
##     aft2$f$fn(aft2$par) # ???

##     ## the events are the same
##     all(event == aft1$f$report(aft1$par)$event) # ok
##     all(event == aft2$f$report(aft2$par)$event) # ok
##     ## The H and logh vectors are very close
##     max(abs(aft1$f$report(aft1$par)$H - aft1$H))       # ok
##     max(abs(aft2$f$report(aft2$par)$H - aft2$H))       # ok
##     max(abs(aft1$f$report(aft1$par)$logh - aft1$logh)) # ok
##     max(abs(aft2$f$report(aft2$par)$logh - aft2$logh)) # ok
##     ##
##     ## the Xs and XDs matrices are very close
##     max(abs(aft1$f$report(aft1$par)$Xs - aft1$Xs))   # ok
##     max(abs(aft1$f$report(aft1$par)$XDs - aft1$XDs)) # ok
##     max(abs(aft2$f$report(aft2$par)$Xs - aft2$Xs))   # ok
##     max(abs(aft2$f$report(aft2$par)$XDs - aft2$XDs)) # ok


##     head(aft1$Xs)
##     head(aft2$Xs)
##     head(aft1$f$report(aft1$par)$Xs)
##     head(aft1$f$report(aft2$par)$Xs)
##     head(aft2$f$report(aft2$par)$Xs)
    
##     aft1$f$gr(aft1$par) # ok
##     delta=function(i,eps=1e-5) {new=rep(0,5); new[i]=eps; new}
##     for (i in 1:length(aft1$par))
##         print((aft1$f$fn(aft1$par+delta(i))-aft1$f$fn(aft1$par-delta(i)))/(2e-5))
##     ## Gradient at the 'final' value are not the same
##     for (i in 1:length(aft2$par))
##         print((aft2$negll(aft1$par+delta(i))-aft2$negll(aft1$par-delta(i)))/(2e-5))

##     ## gradient at the initial value are the same
##     aft1$f$gr(aft1$init) # ok
##     delta=function(i,eps=1e-5) {new=rep(0,5); new[i]=eps; new}
##     for (i in 1:length(aft1$par))
##         print((aft1$f$fn(aft1$init+delta(i))-aft1$f$fn(aft1$init-delta(i)))/(2e-5))

##     ## Objective at the initial values are the same
##     as.numeric(aft1$negll(aft1$init)) - aft1$f$fn(aft1$init) # ok
##     as.numeric(aft2$negll(aft2$init)) - aft2$f$fn(aft2$init) # ok

##     ## objectives at the 'final' values are NOT the same
##     as.numeric(aft1$negll(aft1$init+0.1)) - aft1$f$fn(aft1$init+0.1) # ??
##     as.numeric(aft1$negll(aft1$par)) - aft1$f$fn(aft1$par) # ???

##     undebug(aft1$negll)
##     aft1$negll(aft1$par)

##     aft1$init
##     aft1$f$par
##     aft1$f$fn(aft1$init)
##     aft1$f$fn(aft1$f$par)
##     aft1$f$fn(aft1$par)
##     aft1$f$fn(aft2$par)
    
## }


## ## KL using GausLaguerre Quadrature --numerically unstable
##
## KL <- function(object, true_density = "Weibull",
##                relative = TRUE, symmetric = FALSE,
##                shape = 1, scale = 1, shape2 = 1, mixing_par = 0.5,
##                intercept1 = 1, intercept2 = 1, beta = c(1,1),
##                rate =1, location = 0, n_nodes = 110)
##   {
##
##   GaussLaguerreQuadrature <- function(f, n_nodes) {
##
##     library(orthopolynom)
##
##     roots <- sapply(polyroot(laguerre.polynomials(n_nodes)[[n_nodes+1]]), Re)
##
##     weights <- rep(roots/
##       ((n_nodes+1)^2 * as.function(
##     laguerre.polynomials(n_nodes+1)[[n_nodes+2]])(roots)^2), each =
##     length(object@args$event))
##
##     ## transformation
##
##     g <- function(x) f(x*scale)*rep(exp(x)*scale, each = length(object@args$event))
##     g(roots) %*% weights
##   }
##
##   integrand <- function(x) {
##     true_density =  if(true_density == "Weibull") {
##       rep(dweibull(x, shape, scale), each = length(object@args$event))
##     } else if (true_density == "gamma") { rep(dgamma(x, shape, scale = scale),
##                                               each = length(object@args$event))
##     } else if (true_density == "norm") { dnorm(x, location, scale)
##     } else if (true_density == "logis") { dlogis(x, location, scale)
##     } else if (true_density == "mixture Weibull") {
##         ## densities at each quadrature point are evaluated for all covariate
##         ## patterns and then stacked
##       mixing_par *
##         c(t(do.call(rbind,
##                     lapply(x, dweibull, shape = shape,
##                            scale =
##                              exp(intercept1 + object@args$X %*% beta))))) +
##         (1-mixing_par) *
##         c(t(do.call(rbind,
##                     lapply(x, dweibull, shape = shape2,
##                            scale = exp(intercept2 + object@args$X %*% beta)))))
##
##     } else { true_density = true_density # user supplied density
##     }
##
##     newdata <- data.frame(X = do.call(rbind,
##                                       rep(list(object@args$X),length(x))))
##
##     colnames(newdata) <- as.character(
##       tail(as.list(attr(object@args$lm.obj$terms, "variables")),-2))
##     newdata[[object@args$timeVar]] <- rep(x, each = length(object@args$event))
##
##     model_density <-  predict.aft_mixture(object,
##                                    type = "density",
##                                    newdata = newdata)
##
##
##     model_density[model_density == 0] <- .Machine$double.xmin
##
##     if (symmetric) { true_density*log(true_density/model_density) +
##         model_density*log(model_density/true_density)
##     } else if (relative) { -true_density*log(model_density)
##     } else { true_density*log(true_density/model_density) }
##   }
##
##   GaussLaguerreQuadrature(integrand, n_nodes)/length(object@args$event)
## }

## KL using integrate --fast when the number of unique covariate patterns is small
KL_not_vectorized <- function(object, true_density = "Weibull",
                              relative = TRUE, symmetric = FALSE,
                              shape = c(1,1), scale = c(1,1), rate =1, location = 0,
                              mixing_par = 0.5, beta = c(1,1))
{

    integrand <- function(x, newdata) {
        true_density =  if(true_density == "Weibull") { dweibull(x, shape, scale)
                        } else if (true_density == "gamma") { dgamma(x, shape, scale = scale)
                        } else if (true_density == "norm") { dnorm(x, location, scale)
                        } else if (true_density == "logis") { dlogis(x, location, scale)
                        } else if (true_density == "mixture Weibull") {
                            mixing_par*
                                dweibull(x, shape[1], scale[1]) +
                                (1-mixing_par)*
                                dweibull(x, shape[2], scale[2])
                        } else { true_density = true_density # user supplied density
                        }
        names <- colnames(newdata)
        newdata <- data.frame(X = rep(newdata[[names]],length(x)))
        colnames(newdata) <- as.character(
            tail(as.list(attr(object@args$lm.obj$terms, "variables")),-2))
        newdata[[object@args$timeVar]] <- x
        model_density <-  predict(object,
                                  type = "density",
                                  newdata = newdata)
        if (symmetric) { true_density*log(true_density/model_density) +
                             model_density*log(model_density/true_density)
        } else if (relative) { -true_density*log(model_density)
        } else { true_density*log(true_density/model_density) }
    }
    out <- 0
    ## save unique covariate patterns with duplicate count
    uniqueCov <- data.frame(pre = object@args$X)
    uniqueCov <- aggregate(cbind(uniqueCov[0], nDuplic = 1), uniqueCov, length)
    ## calculate integral only for unique covariate patterns
    for (i in c(1:NROW(uniqueCov))) {
        newdata <- data.frame(X = uniqueCov[i,-NCOL(uniqueCov)])
        colnames(newdata) <-
            as.character(tail(as.list(attr(object@args$lm.obj$terms,"variables")),-2)[[1]])
        out <- out +
            integrate(integrand,
                      0, 30, newdata = newdata)$value*
                                              uniqueCov$nDuplic[i]/sum(uniqueCov$nDuplic)
    }
    return(out)
}
