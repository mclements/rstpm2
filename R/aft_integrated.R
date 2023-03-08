setClass("aft_integrated", representation(args="list"), contains="mle2")

aft_integrated <- function(formula, data, df = 3,
                        tvc = NULL, cure.formula=formula,
                        control = list(parscale = 1, maxit = 1000), init = NULL,
                        weights = NULL, nNodes=20, 
                        timeVar = "", time0Var = "", log.time.transform=TRUE,
                        reltol=1.0e-8, trace = 0, cure = FALSE, mixture = FALSE,
                        contrasts = NULL, subset = NULL, use.gr = TRUE, ...) {
    ## parse the event expression
    eventInstance <- eval(lhs(formula),envir=data)
    stopifnot(length(lhs(formula))>=2)
    eventExpr <- lhs(formula)[[length(lhs(formula))]]
    delayed <- length(lhs(formula))>=4 # indicator for multiple times (cf. strictly delayed)
    surv.type <- attr(eventInstance,"type")
    if (surv.type %in% c("interval","interval2","left","mstate"))
        stop("aft_integrated not implemented for Surv type ",surv.type,".")
    counting <- attr(eventInstance,"type") == "counting"
    ## interval <- attr(eventInstance,"type") == "interval"
    timeExpr <- lhs(formula)[[if (delayed) 3 else 2]] # expression
    if (timeVar == "")
        timeVar <- all.vars(timeExpr)
    ## set up the formulae
    full.formula <- formula
    rhs(full.formula) <- rhs(full.formula) %call+% quote(0)
    if (!is.null(tvc)) {
        tvc.formulas <-
            lapply(names(tvc), function(name)
                call(":",
                     call("as.numeric",as.name(name)),
                     as.call(c(quote(ns),
                               timeExpr,
                               vector2call(list(df=tvc[[name]]))))))
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
    coxph.call$model <- TRUE
    coxph.obj <- eval(coxph.call, envir=parent.frame())
    y <- model.extract(model.frame(coxph.obj),"response")
    data$logHhat <- pmax(-18,log(-log(S0hat(coxph.obj))))
    ## now for the cure fraction
    glm.cure.call = coxph.call
    glm.cure.call[[1]] = as.name("glm")
    glm.cure.call$family = as.name("binomial")
    lhs(glm.cure.call$formula) = as.name("event")
    rhs(glm.cure.call$formula) = rhs(cure.formula)
    ## glm(y ~ X, family=binomial)
    ## browser()
    glm.cure.obj <- eval(glm.cure.call, data)
    Xc = model.matrix(glm.cure.obj, data)
    ##
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
    ## lm0D.obj <- lm(logHhat~nsxD(logtstar,df,intercept=TRUE,cure=cure)-1,dataEvents)
    ## browser()
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
    lpfunc <- function(x,fit,data,var) {
        data[[var]] <- x
        lpmatrix.lm(fit,data)
    }
    ##
    ## surv.type %in% c("right","counting")
    ##
    ## For integrating for time-varying acceleration factors:
    ## - get nodes and weights for Gauss-Legendre quadrature
    ## - get a design matrix for each node
    ## - get a design matrix for the end of follow-up (for the hazard calculations)
    ## - pass that information to C++ for calculation of the integrals and for the hazards
    ## - and we need to do this for the predictions:)
    ##
    gauss = gauss.quad(nNodes)
    ## browser()
    X_list = lapply(1:nNodes, function(i)
        lpmatrix.lm(lm.obj,
                    local({ data[[timeVar]] = (gauss$nodes[i]+1)/2*data[[timeVar]]; data})))
    X <- lpmatrix.lm(lm.obj,data)
    if (delayed && all(time0==0)) delayed <- FALSE # CAREFUL HERE: delayed redefined
    if (delayed) {
        X_list0 = lapply(1:nNodes, function(i)
            lpmatrix.lm(lm.obj,
                        local({ data[[time0Var]] = (gauss$nodes[i]+1)/2*data[[time0Var]]; data})))
    }
    ## Weibull regression
    if (delayed) {
        if (requireNamespace("eha", quietly = TRUE)) {
            survreg1 <- eha::aftreg(formula, data)
            coef1 <- -coef(survreg1) # reversed parameterisation
            coef1 <- coef1[1:(length(coef1)-2)]
        } else coef1 <- rep(0,ncol(X))
    } else {
        survreg1 <- survival::survreg(formula, data)
        coef1 <- coef(survreg1)
        coef1 <- coef1[-1] # assumes intercept included in the formula
    }
    if (ncol(X)>length(coef1)) {
        coef1 <- c(coef1,rep(0,ncol(X) - length(coef1)))
        names(coef1) <- names(coef1b)
    }
    ## browser()
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
    args <- list(init=init,X_list=X_list,
                 X_list0=if (delayed) X_list0 else list(matrix(0,0,0)),
                 wt=wt,event=ifelse(event,1,0),time=time,y=y,
                 time0 = if (delayed) time0 else 0*time,
                 timeVar=timeVar,timeExpr=timeExpr,terms=mt,
                 parscale=parscale, reltol=reltol,
                 Xt=X, Xc= if (mixture) Xc else matrix(0,0,0), maxit=control$maxit,
                 time0=time0, log.time.transform=log.time.transform,
                 trace = as.integer(trace), 
                 boundaryKnots=attr(design,"Boundary.knots"), q.const=t(attr(design,"q.const")),
                 interiorKnots=attr(design,"knots"), design=design, designD=designD,
                 designDD=designDD, cure=as.integer(cure), mixture = as.integer(mixture),
                 data=data, lm.obj = lm.obj, glm.cure.obj = glm.cure.obj, return_type="optim",
                 gweights=gauss$weights, gnodes=gauss$nodes)
    negll <- function(beta) {
        localargs <- args
        localargs$return_type <- "objective"
        localargs$init <- beta
        return(.Call("aft_integrated_model_output", localargs, PACKAGE="rstpm2"))
    }
    gradient <- function(beta) {
        localargs <- args
        localargs$return_type <- "gradient"
        localargs$init <- beta
        return(as.vector(.Call("aft_integrated_model_output", localargs, PACKAGE="rstpm2")))
    }
    parnames(negll) <- names(init)
    args$negll = negll
    args$gradient = gradient
    ## MLE
    if (delayed && use.gr) { # initial search using nmmin (conservative -- is this needed?)
        args$return_type <- "nmmin"
        args$maxit <- 50
        fit <- .Call("aft_integrated_model_output", args, PACKAGE="rstpm2")
        args$maxit <- control$maxit
    }
    optim_step <- function(use.gr) {
        args$return_type <<- if (use.gr) "vmmin" else "nmmin"
        fit <- .Call("aft_integrated_model_output", args, PACKAGE="rstpm2")
        coef <- as.vector(fit$coef)
        hessian <- fit$hessian
        names(coef) <- rownames(hessian) <- colnames(hessian) <- names(init)
        args$init <<- coef
        ## we could use mle2() to calculate vcov by setting eval.only=FALSE
        mle2 <- if (use.gr) bbmle::mle2(negll, coef, vecpar=TRUE, control=control,
                                        gr=gradient, ..., eval.only=TRUE)
                else bbmle::mle2(negll, coef, vecpar=TRUE, control=control, ..., eval.only=TRUE)
        ## browser()
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
    ## browser()
    ## mle2 <-  bbmle::mle2(negll, init, vecpar=TRUE, control=control, ...)
    mle2 <- optim_step(use.gr)
    if (all(is.na(mle2@vcov)) && use.gr) {
        args$init <- init
        mle2 <- optim_step(FALSE)
    }
    out <- as(mle2, "aft_integrated")
    out@args <- args
    attr(out,"nobs") <- length(out@args$event) # for logLik method
    return(out)
}
