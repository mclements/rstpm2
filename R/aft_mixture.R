setClass("aft_mixture", representation(args="list"), contains="mle2")

aft_mixture <- function(formula, data, smooth.formula = NULL, df = 3,
                        tvc = NULL, cure.formula=formula,
                        control = list(parscale = 1, maxit = 1000), init = NULL,
                        weights = NULL, 
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
        stop("aft_mixture not implemented for Surv type ",surv.type,".")
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
                               call("log",timeExpr),
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
    ## glm.cure.obj <- eval(glm.cure.call, envir=parent.frame())
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
    ## initialise values
    ind0 <- FALSE
    map0 <- 0L
    which0 <- 0
    wt0 <- 0
    ## ttype <- 0
    ## surv.type %in% c("right","counting")
    ##
    ## For integrating for time-varying acceleration factors:
    ## - get nodes and weights for Gauss-Legendre quadrature
    ## - get a design matrix for each node
    ## - get a design matrix for the end of follow-up (for the hazard calculations)
    ## - pass that information to C++ for calculation of the integrals and for the hazards
    ## - and we need to do this for the predictions:)
    ##
    X <- lpmatrix.lm(lm.obj,data)
    ## Xc <- model.matrix(coxph.obj, data)
    XD <- grad1(lpfunc,data[[timeVar]],lm.obj,data,timeVar,log.transform=log.time.transform)
    XD <- matrix(XD,nrow=nrow(X))
    Xc0 <- XD0 <- X0 <- matrix(0,1,ncol(X))
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
        Xc0 = lpmatrix.lm(glm.cure.obj, data0)
        wt0 <- wt[ind0]
        XD0 <- grad1(lpfunc,data0[[timeVar]],lm.obj,data0,timeVar,log.transform=log.time.transform)
        XD0 <- matrix(XD0,nrow=nrow(X0))
        rm(data0)
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
        coef1 <- coef1[-1] # assumes intercept included in the formula; ignores smooth.formula
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
    args <- list(init=init,X=X,XD=XD,wt=wt,event=ifelse(event,1,0),time=time,y=y,
                 timeVar=timeVar,timeExpr=timeExpr,terms=mt,
                 delayed=delayed, X0=X0, XD0=XD0, Xc0=Xc0, wt0=wt0, parscale=parscale, reltol=reltol,
                 Xc=if (mixture) Xc else matrix(0,0,0), maxit=control$maxit,
                 time0=if (delayed) time0[time0>0] else NULL, log.time.transform=log.time.transform,
                 trace = as.integer(trace), map0 = map0 - 1L, ind0 = ind0, which0 = which0 - 1L,
                 boundaryKnots=attr(design,"Boundary.knots"), q.const=t(attr(design,"q.const")),
                 interiorKnots=attr(design,"knots"), design=design, designD=designD,
                 designDD=designDD, cure=as.integer(cure), mixture = as.integer(mixture),
                 data=data, lm.obj = lm.obj, glm.cure.obj = glm.cure.obj, return_type="optim")
    negll <- function(beta) {
        localargs <- args
        localargs$return_type <- "objective"
        localargs$init <- beta
        return(.Call("aft_mixture_model_output", localargs, PACKAGE="rstpm2"))
    }
    gradient <- function(beta) {
        localargs <- args
        localargs$return_type <- "gradient"
        localargs$init <- beta
        return(as.vector(.Call("aft_mixture_model_output", localargs, PACKAGE="rstpm2")))
    }
    parnames(negll) <- names(init)
    args$negll = negll
    args$gradient = gradient
    ## MLE
    if (delayed && use.gr) { # initial search using nmmin (conservative -- is this needed?)
        args$return_type <- "nmmin"
        args$maxit <- 50
        fit <- .Call("aft_mixture_model_output", args, PACKAGE="rstpm2")
        args$maxit <- control$maxit
    }
    optim_step <- function(use.gr) {
        args$return_type <<- if (use.gr) "vmmin" else "nmmin"
        fit <- .Call("aft_mixture_model_output", args, PACKAGE="rstpm2")
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
    out <- as(mle2, "aft_mixture")
    out@args <- args
    attr(out,"nobs") <- length(out@args$event) # for logLik method
    return(out)
}

setMethod("nobs", "aft_mixture", function(object, ...) length(object@args$event))

predict.aft_mixture =  function(object,newdata=NULL,
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
        Xc <- lpmatrix.lm(args$glm.cure.obj, newdata)
        time <- eval(args$timeExpr,newdata)
    }
    if (type %in% c("hr","sdiff","hdiff","meansurvdiff","or","af","accfac")) {
        newdata2 <- exposed(newdata)
        X2 <- lpmatrix.lm(args$lm.obj, newdata2)
        XD2 <- grad1(lpfunc,newdata2[[object@args$timeVar]],
                     log.transform=object@args$log.time.transform)
        XD2 <- matrix(XD2,nrow=nrow(X))
        time2 <- eval(args$timeExpr,newdata2) # is this always equal to time?
        Xc2 = model.matrix(args$glm.cure.obj, newdata2)
    }
    if (type == "gradh") {
        return(predict.aft_mixture.ext(object, type="gradh", time=time, X=X, XD=XD))
    }
    ## colMeans <- function(x) colSums(x)/apply(x,2,length)
    local <-  function (object, newdata=NULL, type="surv", exposed, ...)
    {
        args <- object@args
        betafull <- coef(object)
        beta <- betafull[1:ncol(args$X)]
        betac <- betafull[(ncol(args$X)+1):(ncol(args$X)+ncol(args$Xc))]
        betas <- betafull[-(1:(ncol(args$X)+ncol(args$Xc)))]
                                        #beta <- betafull[1:ncol(args$X)]
                                        #betas <- betafull[-(1:ncol(args$X))]
        tt <- args$terms
        eta <- as.vector(X %*% beta)
        logtstar <- log(time) - eta
        etas <- as.vector(predict(args$design, logtstar) %*% betas)
        H <- exp(etas)
        S <- exp(-H)
                                        #browser()
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
            accfac <- eta - log(1-time*etaD)
            accfac2 <- eta2 - log(1-time2*etaD2)
            return(exp(accfac2-accfac))
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

setMethod("predict", "aft_mixture", predict.aft_mixture)
cloglog <- function(x) log(-log(x))
cexpexp <- function(x) exp(-exp(x))

## setMethod("predictnl", "aft_mixture",
##           function(object,fun,newdata=NULL,link=c("I","log","cloglog","logit"), gd=NULL, ...)
##           {
##             link <- match.arg(link)
##             linkf <- eval(parse(text=link))
##             if (is.null(newdata) && !is.null(object@args$data))
##               newdata <- object@args$data
##             localf <- function(coef,...)
##             {
##               object@args$init <- object@fullcoef <- coef
##               linkf(fun(object,...))
##             }
##             numDeltaMethod(object,localf,newdata=newdata,gd=gd,...)
##           })

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
##     ## browser()
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

plot.aft_mixture.meansurv <- function(x, y=NULL, times=NULL, newdata=NULL, type="meansurv", exposed=NULL, add=FALSE, ci=!add, rug=!add, recent=FALSE,
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
plot.aft_mixture.base <-
    function(x,y,newdata=NULL,type="surv",
             xlab=NULL,ylab=NULL,line.col=1,ci.col="grey",lty=par("lty"),
             add=FALSE,ci=!add,rug=!add,
             var=NULL,exposed=incrVar(var),times=NULL,...) {
        if (type %in% c("meansurv","meansurvdiff","af")) {
            return(plot.aft_mixture.meansurv(x,times=times,newdata=newdata,type=type,xlab=xlab,ylab=ylab,line.col=line.col,ci.col=ci.col,
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
setMethod("plot", signature(x="aft_mixture", y="missing"),
          function(x,y,newdata=NULL,type="surv",
                   xlab=NULL,ylab=NULL,line.col=1,ci.col="grey",lty=par("lty"),
                   add=FALSE,ci=!add,rug=!add,
                   var=NULL,exposed=incrVar(var),times=NULL,...)
              plot.aft_mixture.base(x=x, y=y, newdata=newdata, type=type, xlab=xlab,
                                    ylab=ylab, line.col=line.col, lty=lty, add=add,
                                    ci=ci, rug=rug, var=var, exposed=exposed, times=times, ...)
          )
predict.aft_mixture.ext <- function(obj, type=c("survival","haz","gradh"),
                                    time=obj@args$time, X=obj@args$X, XD=obj@args$XD) {
    type <- match.arg(type)
    localargs <- obj@args
    localargs$return_type <- type
    localargs$X <- X
    localargs$XD <- XD
    localargs$time <- time
    as.matrix(.Call("aft_mixture_model_output", localargs, PACKAGE="rstpm2"))
}
