setClass("aft_integrated", representation(args="list"), contains="mle2")

aft_integrated <- function(formula, data, df = 3,
                        tvc = NULL, cure.formula=formula,
                        control = list(parscale = 1, maxit = 1000), init = NULL,
                        weights = NULL, nNodes=20, tvc.intercept=TRUE,
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
    loglpfunc <- function(x,fit,data,var) {
      data[[var]] <- exp(x)
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
                        local({ data[[timeVar]] = (gauss$nodes[i]+1)/2*data[[time0Var]]; data})))
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
    negll2 <- function(beta) {
        localargs <- args
        localargs$return_type <- "objective2"
        localargs$init <- beta
        return(.Call("aft_integrated_model_output", localargs, PACKAGE="rstpm2"))
    }
    gradient2 <- function(beta) {
        localargs <- args
        localargs$return_type <- "gradient2"
        localargs$init <- beta
        return(as.vector(.Call("aft_integrated_model_output", localargs, PACKAGE="rstpm2")))
    }
    parnames(negll2) <- parnames(negll) <- names(init)
    args$negll = negll
    args$gradient = gradient
    args$negll2 = negll2
    args$gradient2 = gradient2
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

setMethod("nobs", "aft_integrated", function(object, ...) length(object@args$event))

setMethod("predict", "aft_integrated",
          function(object,newdata=NULL,
                   type=c("surv","cumhaz","hazard","density","hr","sdiff","hdiff","loghazard","link","meansurv","meansurvdiff","odds","or","meanhaz","af","fail","accfac","gradh"),
                   grid=FALSE,seqLength=300,level=0.95,
                   se.fit=FALSE,link=NULL,exposed=incrVar(var),var=NULL,keep.attributes=TRUE,...) {
              predict.aft_integrated(object, newdata, type, grid, seqLength, level, se.fit, link, exposed, var, keep.attributes, ...)
          })

predict.aft_integrated =  function(object,newdata=NULL,
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
        X <- lpmatrix.lm(args$lm.obj, newdata)
        XD <- grad1(loglpfunc,log(newdata[[object@args$timeVar]]),
                    log.transform=FALSE, newdata=newdata)
        XD <- matrix(XD,nrow=nrow(X))
        Xc <- lpmatrix.lm(args$glm.cure.obj, newdata)
        time <- eval(args$timeExpr,newdata)
        X_list = lapply(args$gnodes, function(gnode)
            lpmatrix.lm(args$lm.obj,
                    local({ newdata[[args$timeVar]] = (gnode+1)/2*newdata[[args$timeVar]]; newdata})))
    }
    if (type %in% c("hr","sdiff","hdiff","meansurvdiff","or","af","accfac")) {
        newdata2 <- exposed(newdata)
        X2 <- lpmatrix.lm(args$lm.obj, newdata2)
        XD2 <- grad1(loglpfunc,log(newdata2[[object@args$timeVar]]),
                     log.transform=FALSE, newdata=newdata2)
        XD2 <- matrix(XD2,nrow=nrow(X))
        time2 <- eval(args$timeExpr,newdata2) # is this always equal to time?
        Xc2 = model.matrix(args$glm.cure.obj, newdata2)
        X_list2 = lapply(args$gnodes, function(gnode)
            lpmatrix.lm(args$lm.obj,
                    local({ newdata2[[args$timeVar]] = (gnode+1)/2*newdata2[[args$timeVar]]; newdata2})))
    }
    if (type == "gradh") {
        return(predict.aft_integrated.ext(object, type="gradh", time=time, X=X, XD=XD,
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
        ## browser()
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
        logtstar <- log(timestar)
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
            return(eta)
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

predict.aft_integrated.ext <- function(obj, type=c("survival","haz","gradh"),
                                       time=obj@args$time, X=obj@args$X, XD=obj@args$XD,
                                       X_list=obj@args$X_list, Xc=obj@args$Xc) {
    type <- match.arg(type)
    localargs <- obj@args
    localargs$return_type <- type
    localargs$X <- X
    localargs$XD <- XD
    localargs$Xc <- Xc
    localargs$X_list <- X_list
    localargs$time <- time
    as.matrix(.Call("aft_integrated_model_output", localargs, PACKAGE="rstpm2"))
}

setMethod("plot", signature(x="aft_integrated", y="missing"),
          function(x,y,newdata=NULL,type="surv",
                   xlab=NULL,ylab=NULL,line.col=1,ci.col="grey",lty=par("lty"),
                   add=FALSE,ci=!add,rug=!add,
                   var=NULL,exposed=incrVar(var),times=NULL,...)
              plot.aft.base(x=x, y=y, newdata=newdata, type=type, xlab=xlab,
                              ylab=ylab, line.col=line.col, lty=lty, add=add,
                              ci=ci, rug=rug, var=var, exposed=exposed, times=times, ...)
          )
