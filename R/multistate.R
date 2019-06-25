markov_msm <- 
    function (x, trans, t = c(0,1), newdata = NULL, init=NULL,
              tmvar = NULL, 
              sing.inf = 1e+10, method="adams", rtol=1e-10, atol=1e-10, slow=FALSE,
              min.tm=1e-8,
              utility=function(t) rep(1, nrow(trans)),
              use.costs=FALSE,
              transition.costs=lapply(1:sum(!is.na(trans)), function(i) function(t) 1000), # per transition
              state.costs=function(t) rep(0,nrow(trans)), # per unit time
              discount.rate = 0,
              block.size=500,
              debug = FALSE,
              ...)
{
    call <- match.call()
    inherits <- function(x, ...)
        base::inherits(x, ...) ||
            (base::inherits(x, c("zeroModel","hrModel","stratifiedModel"))
                && base::inherits(x$base, ...))
    base.classes <- c("stpm2","pstpm2","glm","survPen","gam","aft","survreg")
    stopifnot(all(sapply(x, function(xi) inherits(xi,base.classes) | is.function(xi))))
    stopifnot(!is.null(newdata))
    stopifnot(sum(!is.na(trans)) == length(x))
    if (is.null(init)) init <- "[<-"(rep(0,nrow(trans)),1,1)
    stopifnot(length(init) == nrow(trans))
    ## if newdata are many, then separate into blocks
    if (nrow(newdata)>block.size) {
        n.blocks <- ceiling(nrow(newdata)/block.size)
        indices <- parallel::splitIndices(nrow(newdata), n.blocks)
        lst <-
            parallel::mclapply(indices,
                               function(index)
                                   markov_msm(x = x, trans = trans, t = t,
                                              newdata = newdata[index,],
                                              init = init, tmvar = tmvar,
                                              sing.inf = sing.inf,
                                              method = method, rtol = rtol,
                                              atol = atol, slow = slow,
                                              min.tm = min.tm, utility = utility,
                                              use.costs = use.costs,
                                              transition.costs = transition.costs,
                                              state.costs = state.costs,
                                              discount.rate = discount.rate,
                                              block.size = block.size,
                                              ...))
        return(do.call(rbind.markov_msm,lst))
    }
    if (use.costs)
        slow <- TRUE
    nt <- length(t)
    if (nt < 2) 
        stop("number of times should be at least two")
    stopifnot(length(utility(t[2])) %in% c(1,nrow(trans)))
    if (is.null(tmvar) && all(sapply(x,inherits,c("stpm2","pstpm2","aft","survPen"))))
        tmvar <- sapply(x,function(object)
            if(inherits(object,c("stpm2","pstpm2"))) object@timeVar
            else if (inherits(object,"aft")) object@args$timeVar
            else object$t1.name)
    stopifnot(!is.null(tmvar))
    stopifnot(length(tmvar) %in% c(1,length(x)))
    if (length(tmvar)==1) tmvar <- rep(tmvar, length(x))
    for(i in 1:length(x)) {
        if(inherits(x[[i]],"survPen") && !inherits(x[[i]],"survPenWrap"))
            x[[i]] <- survPenWrap(x[[i]])
        if(inherits(x[[i]],"gam") && !inherits(x[[i]],"gamWrap"))
            x[[i]] <- gamWrap(x[[i]])
        else if(is.function(x[[i]]))
            x[[i]] <- hazFun(x[[i]], tmvar[i])
    }
    ntr <- sum(!is.na(trans))
    nobs <- nrow(newdata)
    nstates <- nrow(trans)
    n <- nstates*nobs
    cs <- sapply(x,function(item) length(coef(item)))
    cumcs <- cumsum(c(0,cs))
    ncoef <- sum(cs)
    lambda.discount <- log(1+discount.rate)
    matwhich <- function(ind) {
        stopifnot(sum(ind,na.rm=TRUE) == 1)
        nr <- nrow(ind)
        nc <- ncol(ind)
        wind <- which(ind)-1
        c(wind %% nr + 1, wind %/% nr + 1)
    }
    coef.surv_list <- function(objects)
        sapply(objects,function(x) coef(x)) # sapply(objects, coef) FAILS
    vcov.surv_list <- function(objects) {
        vcovs <- lapply(objects, function(x) vcov(x))
        lengths <- sapply(vcovs,nrow)
        index <- cumsum(c(0,lengths))
        m <- matrix(0,sum(lengths),sum(lengths))
        for (i in 1:length(objects))
            m[(index[i]+1):index[i+1],(index[i]+1):index[i+1]] <- vcovs[[i]]
        rownames(m) <- colnames(m) <- names(coef(objects))
        m
    }
    dp <- function(t, y, parms, ...) {
        P <- matrix(y[1:(nstates*nobs)],nstates,nobs) # only one state vector per individual
        Pu <- array(y[nstates*nobs+1:(nstates*ncoef*nobs)],c(nstates,ncoef,nobs))
        QC <- Q <- array(0, c(nstates, nstates, nobs))
        QCu <- Qu <- array(0, c(nstates, nstates, ncoef, nobs))
        for (i in 1:ntr) {
            newdata[[tmvar[i]]] <- pmax(min.tm,t) # special case: t=0
            ij <- matwhich(trans == i)
            from <- ij[1]
            to <- ij[2]
            Q[from,to,] <- predict(x[[i]], newdata=newdata, type="haz")
            Qu[from,to,1:cs[i]+cumcs[i],] <- t(predict(x[[i]], newdata=newdata, type="gradh")) # transposed
            if (use.costs) {
                QC[from,to,] <- Q[from,to,]*transition.costs[[i]](t)
                QCu[from,to,1:cs[i]+cumcs[i],] <- Qu[from,to,1:cs[i]+cumcs[i],]*transition.costs[[i]](t)
            }
        }
        Q[is.na(Q)] <- Qu[is.na(Qu)] <- QC[is.na(QC)] <- QCu[is.na(QCu)] <- 0
        Q[is.infinite(Q) & Q > 0] <- Qu[is.infinite(Qu) & Qu > 0] <- sing.inf
        QC[is.infinite(QC) & QC > 0] <- QCu[is.infinite(QCu) & QCu > 0] <- sing.inf
        if (slow) {
            dTCdt <- dPdt <- array(0,c(nstates,nobs))
            dTCudt <- dPudt <- array(0,c(nstates,ncoef,nobs))
            for (i in 1:nobs) {
                Qi <- Q[,,i]                  # nstates*nstates
                Pi <- P[,i]                   # nstates
                diag(Qi) <- -rowSums(Qi)
                dPdt[,i] <- t(Pi)%*%Qi        # nstates
                Qui <- Qu[,,,i]               # nstates*nstates*ncoef
                Pui <- Pu[,,i]                # nstates*ncoef
                for (j in 1:ncoef)
                    diag(Qui[,,j]) <- -rowSums(Qui[,,j])
                dPudt[,,i] <- t(t(Pui)%*%Qi) + apply(Qui,3,function(x) t(Pi)%*%x) # nstates*ncoef
                if (use.costs) {
                    QCi <- QC[,,i]                # nstates*nstates
                    QCui <- QCu[,,,i]             # nstates*nstates*ncoef
                    dTCdt[,i] <- t(Pi)%*%QCi + Pi*state.costs(t)      # nstates
                    dTCudt[,,i] <- t(t(Pui)%*%QCi) + apply(QCui,3,function(x) t(Pi)%*%x) # nstates*ncoef
                }
            }
        } else {
            Qu2 <- array(Qu,c(nstates*nstates,ncoef,nobs))
            QCu2 <- array(QCu,c(nstates*nstates,ncoef,nobs))
            lst <- .Call("multistate_ddt",P,Pu,Q,Qu2,PACKAGE="rstpm2")
            dPdt <- lst[[1]]
            dPudt <- lst[[2]]
            dTCdt <- array(0,c(nstates,nobs))        # not implemented!
            dTCudt <- array(0,c(nstates,ncoef,nobs)) # not implemented!
        }
        if (use.costs)
            list(c(dPdt,
                   dPudt,
                   P*utility(t)*exp(-lambda.discount*t),
                   Pu*utility(t)*exp(-lambda.discount*t),
                   dTCdt*exp(-lambda.discount*t),
                   dTCudt*exp(-lambda.discount*t)
                   ))
        else
            list(c(dPdt,
                   dPudt,
                   P*utility(t)*exp(-lambda.discount*t),
                   Pu*utility(t)*exp(-lambda.discount*t)
                   ))
    }
    init <- if (use.costs) c(rep(init,nobs),rep(0,2*nobs*nstates+3*nobs*nstates*ncoef))
            else
                c(rep(init,nobs),rep(0,nobs*nstates+2*nobs*nstates*ncoef))
    res <- ode(y = init, times = t, func = dp, method=method, atol=atol, rtol=rtol,
               ...)
    class(x) <- "surv_list"
    Sigma <- vcov(x)
    P <- array(res[,1+1:(nstates*nobs)],c(nt,nstates,nobs))
    Pu <- array(res[,1+nobs*nstates+1:(nobs*nstates*ncoef)],c(nt,nstates,ncoef,nobs))
    L <- array(res[,1+nstates*nobs+nstates*nobs*ncoef+1:(nstates*nobs)],c(nt,nstates,nobs))
    Lu <- array(res[,1+2*nstates*nobs+nstates*nobs*ncoef+1:(nobs*nstates*ncoef)],c(nt,nstates,ncoef,nobs))
    state.names <- if(is.null(rownames(trans))) 1:nrow(trans) else rownames(trans)
    dimnames(P) <- dimnames(L) <- 
        list(time=t, state= state.names, obs=rownames(newdata))
    coef.names <- unlist(lapply(x, function(xi) names(coef(xi))))
    dimnames(Pu) <- dimnames(Lu) <-
        list(time=t,
             state=state.names,
             coef=coef.names,
             obs=rownames(newdata))
    structure(list(time = t, P=P, Pu=Pu, L=L, Lu=Lu, 
                   res=if (debug) res else NULL, # for debugging only
                   vcov=Sigma, newdata=newdata, trans=trans,
                   call=call),
              class="markov_msm")
}

hazFun <- function(f, tmvar="t", ...) {
    if(all(c("newdata","t") %in% names(formals(f))))
        newf <- function(newdata) f(t=newdata[[tmvar]],newdata=newdata,...)
    else if("t" %in% names(formals(f)))
        newf <- function(newdata) f(t=newdata[[tmvar]], ...)
    else if("newdata" %in% names(formals(f)))
        newf <- function(newdata) f(newdata=newdata, ...)
    else 
        newf <- function(newdata) f(...)
    out <- list(haz=newf)
    class(out) <- "hazFun"
    out
}
predict.hazFun <- function(object, newdata, type=c("haz","gradh")) {
    type <- match.arg(type)
    val <- if (type=="haz") object$haz(newdata) else 0
    if (length(val)==1 && length(val)<nrow(newdata))
        val <- rep(val,nrow(newdata))
    val
}
coef.hazFun <- function(object, ...) c(hazFun=0)
vcov.hazFun <- function(object, ...) matrix(0,1,1,FALSE,list("hazFun","hazFun"))

as.data.frame.markov_msm <- function(x, row.names=NULL, optional=FALSE,
                                     ci=TRUE,
                                     P.conf.type="log-log", L.conf.type="log",
                                     P.range=c(0,1), L.range=c(0,Inf),
                                     ...) {
    state.names <- rownames(x$trans)
    perm <- function(x) aperm(x, c(2,3,1))
    seP <- apply(x$Pu, c(1,2), function(m) sqrt(colSums(m * (vcov(x) %*% m)))) # c(nobs,nt,nstates)
    seP <- array(seP, dim(x$Pu)[c(4,1,2)]) # c(nobs,nt,nstates))
    seP <- aperm(seP,c(2,3,1))
    seL <- apply(x$Lu, c(1,2), function(m) sqrt(colSums(m * (vcov(x) %*% m)))) # c(nobs,nt,nstates)
    seL <- array(seL, dim(x$Lu)[c(4,1,2)]) # c(nobs,nt,nstates))
    seL <- aperm(seL,c(2,3,1))
    out <- data.frame(P=as.vector(perm(x$P)),
                      seP=as.vector(perm(seP)),
                      L=as.vector(perm(x$L)),
                      seL=as.vector(perm(seL)))
    dimnames. <- dimnames(perm(x$P)) 
    dimnames.$obs <- 1:length(dimnames.$obs)
    y <- expand.grid(dimnames.)
    y$time <- as.numeric(attr(y$time,"levels"))[y$time]
    out <- cbind(y,out)
    if (ci) {
        tmp <- surv.confint(out$P,out$seP, conf.type=P.conf.type, min.value=P.range[1], max.value=P.range[2])
        out$P.lower=tmp$lower
        out$P.upper=tmp$upper
        tmp <- surv.confint(out$L,out$seL, conf.type=L.conf.type, min.value=L.range[1], max.value=L.range[2])
        out$L.lower=tmp$lower
        out$L.upper=tmp$upper
    }
    newdata <- as.data.frame(x$newdata[as.numeric(out$obs), , drop=FALSE])
    names(newdata) <- names(x$newdata)
    out <- "rownames<-"(cbind(newdata, out),rownames(out))
    if(!is.null(row.names)) rownames(out) <- row.names
    out
}
coef.markov_msm <- function(object, ...)
    sapply(object,coef)
standardise <- function(x, ...)
    UseMethod("standardise")
standardise.markov_msm <- function(x,
                                   weights = rep(1,nrow(x$newdata)),
                                   normalise = TRUE, ...) {
    call <- match.call
    state.names <- rownames(x$trans)
    add.dim <- function(x) array(x,dim=c(dim(x),1))
    if (length(weights)==0) weights <- 1
    stopifnot(length(weights) %in% c(1,nrow(x$newdata)))
    if (normalise) weights <- weights/sum(weights)
    wtd.mean <- function(x) sum(x*weights)
    Sigma <- vcov(x)
    P <- add.dim(apply(x$P,c(1,2),wtd.mean))
    Pu <- add.dim(apply(x$Pu,1:3,wtd.mean))
    L <- add.dim(apply(x$L,c(1,2),wtd.mean))
    Lu <- add.dim(apply(x$Lu,1:3,wtd.mean))
    dimnames(P) <- dimnames(L) <- list(time=x$time, state=state.names, obs=1)
    oldnames <- dimnames(x$Pu)
    dimnames(Pu) <- dimnames(Lu) <- list(time=x$time, state=state.names, coef=oldnames[[3]], obs=1)
    df <- x$newdata[0,]
    out <- list(P=P, Pu=Pu, L=L, Lu=Lu, vcov=vcov(x), time=x$time, newdata=df, trans=x$trans,
                call=call)
    class(out) <- "markov_msm"
    out
}
vcov.survPen <- function(object) object$Vp
"coef<-.survPen" <- function(this,value) {
    this$coefficients <- value
    this
}
survPenWrap <- function(object) {
    stopifnot(inherits(object, "survPen"))
    if (!inherits(object, "survPenWrap"))
        class(object) <- c("survPenWrap", class(object))
    object
}
predict.survPenWrap <- function(object, newdata, type=NULL, ...) {
    if (is.null(type)) NextMethod("predict", object)
    else if (type=="haz") predict(object, newdata=newdata, do.surv=FALSE)$haz
    else if (type=="gradh") predict(object, newdata=newdata, type="lpmatrix") *
                                predict(object, newdata=newdata, do.surv=FALSE)$haz
    else NextMethod("predict", object)
}
predict.glm <- function (object, newdata=NULL, type=NULL, ...) {
    if (is.null(type)) NextMethod("predict", object)
    if(type=="gradh" && object$family$link != "log")
        stop("Currently only implemented for a log link")
    ## stopifnot() # Poisson family with log link?
    ## assumes response is a rate
    if(type=="haz") stats::predict.glm(object, newdata=newdata, type="response") 
    else if (type=="gradh")
        stats::predict.glm(object, newdata=newdata, type="response") *
            lpmatrix.lm(object, newdata=newdata)
    else NextMethod("predict", object)
}
gamWrap <-  function(object) {
    stopifnot(inherits(object, "gam"))
    if (!inherits(object, "gamWrap"))
        class(object) <- c("gamWrap", class(object))
    object
}
predict.gamWrap <- function (object, newdata=NULL, type=NULL, ...) {
    if (is.null(type)) NextMethod("predict", object)
    if(type=="gradh" && object$family$link != "log")
        stop("Currently only implemented for a log link")
    ## stopifnot() # Poisson family with log link?
    ## assumes response is a rate
    if (type=="haz") predict(object, newdata=newdata, type="response") 
    else if (type=="gradh") {
        as.vector(predict(object, newdata=newdata, type="response")) *
            predict(object, newdata=newdata, type="lpmatrix")
    } else NextMethod("predict", object)
}
surv.confint <- function(p, se, conf.type=c("log","log-log","plain","logit","arcsin"),
                         conf.int=0.95, min.value=0, max.value=1) {
    conf.type <- match.arg(conf.type)
    zval <- qnorm(1 - (1 - conf.int)/2)
    selog <- se/p # transform to selog
    if (conf.type == "plain") {
        se2 <- selog * p * zval
        list(lower = pmax(p - se2, min.value), upper = pmin(p + se2, max.value))
    }
    else if (conf.type == "log") {
        xx <- ifelse(p == 0, NA, p)
        selog2 <- zval * selog
        temp1 <- exp(log(xx) - selog2)
        temp2 <- exp(log(xx) + selog2)
        list(lower = temp1, upper = pmin(temp2, max.value))
    }
    else if (conf.type == "log-log") {
        xx <- ifelse(p == 0 | p == 1, NA, p)
        selog2 <- zval * selog/log(xx)
        temp1 <- exp(-exp(log(-log(xx)) - selog2))
        temp2 <- exp(-exp(log(-log(xx)) + selog2))
        list(lower = temp1, upper = temp2)
    }
    else if (conf.type == "logit") {
        xx <- ifelse(p == 0, NA, p)
        selog2 <- zval * selog * (1 + xx/(1 - xx))
        temp1 <- 1 - 1/(1 + exp(log(p/(1 - p)) - selog2))
        temp2 <- 1 - 1/(1 + exp(log(p/(1 - p)) + selog2))
        list(lower = temp1, upper = temp2)
    }
    else if (conf.type == "arcsin") {
        xx <- ifelse(p == 0, NA, p)
        selog2 <- 0.5 * zval * selog * sqrt(xx/(1 - xx))
        list(lower = (sin(pmax(0, asin(sqrt(xx)) - selog2)))^2, 
             upper = (sin(pmin(pi/2, asin(sqrt(xx)) + selog2)))^2)
    }
}
print.markov_msm <- function(x, 
                             digits=5,
                             se=FALSE, ci=TRUE,
                             P.conf.type="log-log", L.conf.type="log",
                             P.range=c(0,1), L.range=c(0,Inf),
                             ...) {
    df <- as.data.frame(x, ci=ci,
                        P.conf.type=P.conf.type, L.conf.type=L.conf.type,
                        P.range=P.range, L.range=L.range)
    if (!se) df$seP <- df$seL <- NULL
    print(df, digits=digits, ...)
    invisible(x)
}

zeroModel <- function(object) {
    if (inherits(object, "zeroModel")) object
    else structure(list(base=object), class="zeroModel")
}
predict.zeroModel <- function(object, type=c("haz","gradh"), ...) {
    type <- match.arg(type)
    0.0*predict(object$base, type=type...) 
}
## we keep the same length/dim to support comparisons between different interventions
coef.zeroModel <- function(object, ...) 0.0*coef(object$base, ...)
vcov.zeroModel <- function(object, ...) 0.0*vcov(object$base, ...)

hrModel <- function(object, hr=1, ci=NULL) {
    stopifnot(is.null(ci) || (is.numeric(ci) && length(ci)==2))
    stopifnot(is.numeric(hr) && length(hr)==1 && hr>0)
    newobject <- if (inherits(object, "hrModel")) object else structure(list(base=object),
                                                                        class="hrModel")
    attr(newobject, "loghr") <- log(hr)
    attr(newobject, "seloghr") <- if (is.null(ci)) 0 else log(ci[2]/ci[1])/2/qnorm(0.975)
    newobject
}
predict.hrModel <- function(object, type=c("haz","gradh"), ...) {
    type <- match.arg(type)
    hr <- exp(attr(object,"loghr"))
    pred1 <- predict(object$base, type="haz", ...)
    if (type=="haz") hr*pred1
    else
        cbind(predict(object$base, type="gradh", ...)*hr, pred1*hr*log(hr))
}
## This has different lengths/dimensions to the base model:
##   wrap base intervention in hrModel(..., hr=1)
coef.hrModel <- function(object, ...) c(coef(object$base), attr(object,"loghr"))
vcov.hrModel <- function(object, ...) {
    out <- rbind(cbind(vcov(object$base),0),0)
    out[nrow(out),ncol(out)] <- attr(object, "seloghr")^2
    out
}

subset.markov_msm <- function(x, subset, ...) {
    e <- substitute(subset)
    r <- eval(e, x$newdata, parent.frame())
    if (!is.logical(r)) 
        stop("'subset' must be logical")
    r <- r & !is.na(r)
    newx <- list(time=x$time,
                 P=x$P[,,r,drop=FALSE],
                 Pu=x$Pu[,,,r,drop=FALSE],
                 ## seP=x$seP[,,r,drop=FALSE],
                 L=x$L[,,r,drop=FALSE],
                 Lu=x$Lu[,,,r,drop=FALSE],
                 ## seL=x$seL[,,r,drop=FALSE],
                 res=NULL, # not strictly needed
                 vcov=x$vcov,
                 newdata=x$newdata[r,,drop=FALSE],
                 trans=x$trans)
    newx$call <- match.call()
    class(newx) <- "markov_msm"
    newx
}
diff.markov_msm <- function(x, y, ...) {
    stopifnot(inherits(x,"markov_msm"))
    stopifnot(inherits(y,"markov_msm"))
    ## vcov, time, trans should all be the same
    stopifnot(all(x$vcov == y$vcov))
    stopifnot(all(x$time == y$time))
    stopifnot(all(x$trans == y$trans, na.rm=TRUE))
    z <- list(P=x$P-y$P, Pu=x$Pu-y$Pu, L=x$L-y$L, Lu=x$Lu-y$Lu)
    z <- c(list(time=x$time,
                vcov=x$vcov,
                trans=x$trans),
           z)
    z$call <- match.call()
    z$newdata <- x$newdata
    z$newdata[x$newdata != y$newdata] <-  NA
    class(z) <- c("markov_msm_diff","markov_msm") # not strictly "markov_msm"...
    z
}
as.data.frame.markov_msm_diff <- function(x, ...)
    as.data.frame.markov_msm(x, ...,
                             P.conf.type="plain", L.conf.type="plain",
                             P.range=c(-Inf, Inf), L.range=c(-Inf, Inf))


plot.markov_msm <- function(x, y, stacked=TRUE, which=c("P","L"),
                            xlab="Time", ylab=NULL, col=2:6, border=col,
                            ggplot2=FALSE, lattice=FALSE, alpha=0.2,
                            strata=NULL,
                            ...) {
    stopifnot(inherits(x,"markov_msm"))
    which <- match.arg(which)
    ## ylab defaults
    if(is.null(ylab)) {
        ylab <- if(which=='P') "Probability" else "Length of stay"
        if(inherits(x,"markov_msm_diff"))
            ylab <- if(which=='P') "Difference in probabilities"
                    else "Difference in lengths of stay"
        if(inherits(x,"markov_msm_ratio"))
            ylab <- if(which=='P') "Ratio of probabilities"
                    else "Ratio of lengths of stay"
    }
    if (ggplot2)
        ggplot.markov_msm(x, which=which, stacked=stacked, xlab=xlab, ylab=ylab,
                          alpha=alpha, ...)
    else if (lattice)
        xyplot.markov_msm(x, which=which, stacked=stacked, xlab=xlab, ylab=ylab,
                          col=col, border=border, strata=strata, ...)
    else {
        if (nrow(x$newdata)>1) {
            warning("More than one set of covariates; covariates have been standardised")
            x <- standardise(x)
        }
        if (!missing(y)) warning("y argument is ignored")
        df <- as.data.frame(x)
        states <- unique(df$state)
        if (stacked) {
            out <- graphics::plot(range(x$time),0:1, type="n", xlab=xlab, ylab=ylab, ...)
            lower <- 0
            for (i in length(states):1) { # put the last state at the bottom
                df2 <- df[df$state==states[i],]
                if (length(lower)==1) lower <- rep(0,nrow(df2))
                upper <- lower+df2$P
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

ggplot.markov_msm <- function(data, mapping=NULL,
                              which=c("P","L"), 
                              stacked = TRUE, alpha=0.2,
                              xlab=NULL, ylab=NULL,
                              ..., environment=parent.frame()) {
    if (requireNamespace("ggplot2", quietly=TRUE)) {
        which <- match.arg(which)
        df <- as.data.frame(data, ci=!stacked)
        if (stacked)
            ggplot2::ggplot(df, if(is.null(mapping))
                                    ggplot2::aes_string(x='time', y=which, fill='state')
                                else mapping) +
                ggplot2::geom_area() + 
                ggplot2::xlab(xlab) + ggplot2::ylab(ylab)
        else {
            lower <- paste0(which,'.lower')
            upper <- paste0(which,'.upper')
            ggplot2::ggplot(df, if (is.null(mapping))
                                     ggplot2::aes_string(x='time', y=which, ymin=lower,
                                                         ymax=upper)
                                 else mapping) +
                 ggplot2::geom_line() + ggplot2::geom_ribbon(alpha=alpha) +
                 ggplot2::facet_grid(stats::reformulate(".","state")) +
                ggplot2::xlab(xlab) + ggplot2::ylab(ylab)
        }
    }
    else stop("Suggested package 'ggplot2' is not available")
}

xyplot.markov_msm <- function(x, data, strata=NULL,
                              which=c("P","L"), 
                              stacked = TRUE, 
                              col=2:6, border=col,
                              ..., environment=parent.frame()) {
    if (requireNamespace("lattice", quietly=TRUE)) {
        which <- match.arg(which)
        ## if (!missing(data) && is.null(strata)) strata <- data # copy from data to strata
        df <- as.data.frame(x, ci=!stacked)
        if (!is.factor(df$state)) df$state <- as.factor(df$state)
        states <- levels(df$state)
        if (stacked) {
            rhs_string <- if (!is.null(strata)) paste('time',deparse(rhs(strata)),sep='|')
                          else 'time'
            lattice::xyplot(stats::reformulate(rhs_string,which), data=df,
                            panel = function(x, y, subscripts, group, ...) {
                                df2 <- df[subscripts,]
                                lattice::panel.xyplot(x,y, type="n")
                                lower <- 0
                                for (i in length(states):1) {
                                    df3 <- df2[df2$state==states[i], , drop=FALSE]
                                    if (length(lower)==1) lower <- rep(0,nrow(df3))
                                    upper <- lower+df3[[which]]
                                    lattice::panel.polygon(c(df3$time,rev(df3$time)),
                                                           c(lower,rev(upper)),
                                                  border=border[i], col=col[i])
                                    lower <- upper
                                }
                            },
                            subscripts=TRUE, ...)
        }
        else {
            lower <- paste0(which,'.lower')
            upper <- paste0(which,'.upper')
            rhs_string <- if (!is.null(strata)) paste('time|state',deparse(rhs(strata)),sep='+')
                          else 'time|state'
            lattice::xyplot(stats::reformulate(rhs_string,which), data=df,
                            panel = function(x, y, subscripts, group, ...) {
                                lattice::panel.xyplot(x,y, type="n")
                                lattice::panel.polygon(c(x,rev(x)),
                                                       c(df[[lower]][subscripts],
                                                         rev(df[[upper]][subscripts])),
                                                       border=NULL, col="grey")
                                lattice::panel.lines(x,y)
                            },
                            subscripts=TRUE, ...)
        }
    }
    else stop("Suggested package 'lattice' is not available")
}

## g(beta) = log(losA) - log(losB)
## g'(beta) = losA'/losA - losB'/losB # which is not defined when los=0
## NB: data are on a log scale!
ratio_markov_msm <- function(x, y, ...) {
    stopifnot(inherits(x,"markov_msm"))
    stopifnot(inherits(y,"markov_msm"))
    ## vcov, time and trans should all be the same
    stopifnot(all(x$vcov == y$vcov))
    stopifnot(all(x$time == y$time))
    stopifnot(all(x$trans == y$trans, na.rm=TRUE))
    z <- list(P=log(x$P/y$P), L=log(x$L/y$L)) # logs of the ratios
    z$Pu <- apply(x$Pu, 3, function(slice) slice/x$P) - apply(y$Pu, 3, function(slice) slice/y$P)
    dim(z$Pu) <- dim(x$Pu)
    dimnames(z$Pu) <- dimnames(x$Pu)
    z$Lu <- apply(x$Lu, 3, function(slice) slice/x$L) - apply(y$Lu, 3, function(slice) slice/y$L)
    dim(z$Lu) <- dim(x$Lu)
    dimnames(z$Lu) <- dimnames(x$Lu)
    z <- c(list(time=x$time,
                vcov=x$vcov,
                trans=x$trans,
                res=NULL),
           z)
    z$call <- match.call()
    z$newdata <- x$newdata
    z$newdata[x$newdata != y$newdata] <-  NA
    class(z) <- c("markov_msm_ratio","markov_msm") # not strictly "markov_msm"...
    z
}
as.data.frame.markov_msm_ratio <- function(x, ...) {
    ## data are on a log scale!!
    z <- as.data.frame.markov_msm_diff(x, ...)
    for (name in c("P","L","P.lower","P.upper","L.lower","L.upper"))
        z[[name]] <- exp(z[[name]])
    z
}

## prev_markov_msm <- function(x, w, ...) {
##     stopifnot(inherits(x,"markov_msm"))
##     stopifnot(nrow(x$trans) == length(w))
##     stopifnot(all(w %in% 0:1))
##     x <- msm2
##     w <- c(1,1,0)
##     pi <- piL <- P <- L <- array(0,dim(x$P))
##     sumP <- sumL <- array(0,dim(x$P)[-2])
##     for (i in 1:length(w)) {
##         P[,i,] <- x$P[,i,]*w[i]
##         L[,i,] <- x$L[,i,]*w[i]
##         sumP <- sumP + P[,i,]
##         sumL <- sumL + L[,i,]
##     }
##     for (i in 1:length(w)) {
##         pi[,i,] <- P[,i,]/sumP
##         piL[,i,] <- P[,i,]/sumL
##     }
##     dw <- diag(w)
##     logit <- stats::poisson()$linkinv
##     z <- list(P=logit(P), L=log(L)) # logit and log
##     z$Pu <- apply(x$Pu, 3, function(slice) slice/x$P) - apply(y$Pu, 3, function(slice) slice/y$P)
##     dim(z$Pu) <- dim(x$Pu)
##     dimnames(z$Pu) <- dimnames(x$Pu)
##     z$Lu <- apply(x$Lu, 3, function(slice) slice/x$L) - apply(y$Lu, 3, function(slice) slice/y$L)
##     dim(z$Lu) <- dim(x$Lu)
##     dimnames(z$Lu) <- dimnames(x$Lu)
##     z <- c(list(time=x$time,
##                 vcov=x$vcov,
##                 trans=x$trans,
##                 res=NULL),
##            z)
##     z$call <- match.call()
##     z$newdata <- x$newdata
##     z$newdata[x$newdata != y$newdata] <-  NA
##     class(z) <- c("markov_msm_prev","markov_msm") # not strictly "markov_msm"...
##     z
## }
## as.data.frame.markov_msm_prev <- function(x, ...) {
##     ## data are on a log scale!!
##     z <- as.data.frame.markov_msm_diff(x, ...)
##     for (name in c("P","L","P.lower","P.upper","L.lower","L.upper"))
##         z[[name]] <- exp(z[[name]])
##     z
## }


bindlast <- function(...) { # bind on last slice for a bag of arrays
    x <- list(...)
    stopifnot(all(sapply(x,is.array)))
    dim1 <- dim(x[[1]])
    stopifnot(all(sapply(x[-1],function(xi) length(dim(xi))==length(dim1))))
    for (i in 1:(length(dim1)-1))
        stopifnot(all(sapply(x[-1],function(xi) dim(xi)[i]==dim1[i])))
    y <- do.call("c",lapply(x,as.vector))
    dims <- dim1
    dimnames. <- dimnames(x[[1]])
    dims[length(dims)] <- sum(sapply(x,function(xi) dim(xi)[length(dims)]))
    if (!is.null(dimnames.[[length(dims)]]))
        dimnames.[[length(dims)]] <- unlist(lapply(x,function(xi) dimnames(xi)[length(dims)]))
    dim(y) <- dims
    dimnames(y) <- dimnames.
    y
}
rbind.markov_msm <- function(..., deparse.level = 1) {
    x <- list(...)
    stopifnot(all(sapply(x,inherits,"markov_msm")))
    ## class, vcov, time and trans should all be the same
    stopifnot(all(sapply(x[-1],function(xi) all(class(xi) == class(x[[1]])))))
    stopifnot(all(sapply(x[-1],function(xi) all(xi$vcov==x[[1]]$vcov))))
    stopifnot(all(sapply(x[-1],function(xi) all(xi$time==x[[1]]$time))))
    stopifnot(all(sapply(x[-1],function(xi) all(xi$trans==x[[1]]$trans, na.rm=TRUE))))
    bind <- function(name) do.call(bindlast, lapply(x,"[[",name))
    newx <- list(time=x[[1]]$time,
                 P=bind("P"),
                 Pu=bind("Pu"),
                 L=bind("L"),
                 Lu=bind("Lu"),
                 vcov=x[[1]]$vcov,
                 newdata=do.call(rbind,lapply(x,"[[", "newdata")),
                 trans=x[[1]]$trans,
                 res=NULL)
    newx$call <- match.call()
    ## is this a good idea? It will be lost if done more than once...
    newx$newdata$.index <- unlist(lapply(1:length(x),
                                         function(i) rep(i,nrow(x[[i]]$newdata))))
    state.names <- rownames(newx$trans)
    class(newx) <- class(x[[1]])
    newx
}
transform.markov_msm <- function(`_data`, ...) {
    `_data`$newdata <- transform(`_data`$newdata, ...)
    `_data`
}
reorder <- function(x,order) {
    l <- attr(x,"levels")
    factor(l[x],l[order])
}
nrow.markov_msm <- function(x, ...) nrow(as.data.frame(x, ...)) 
vcov.markov_msm <- function(object, ...) object$vcov

collapse_markov_msm <- function(object, which=NULL, sep="; ") {
    ## example: which=c(1,2) => combine states 1 and 2
    stopifnot(inherits(object, "markov_msm"))
    stopifnot(!is.null(which))
    if(is.character(which))
        which <- pmatch(which, rownames(object$trans))
    stopifnot(all(which %in% 1:nrow(object$trans)))
    stopifnot(nrow(object$trans)>=2)
    call <- match.call()
    which <- sort(unique(which)) # in case of duplicates
    if(length(which)==1) return(object) # no change
    ## algorithm:
    n.newstates <- nrow(object$trans) - length(which) + 1
    ## initalise
    newstates <- vector("list", n.newstates)
    j <- 0
    base <- NULL
    for (i in 1:nrow(object$trans)) {
        if (i %in% which) {
            if (is.null(base)) {
                j <- j + 1
                base <- j
            }
            newstates[[base]] <- c(newstates[[base]],i)
        }
        else {
            j <- j + 1
            newstates[[j]] <- i
        }
    }
    index <- vector('numeric',nrow(object$trans))
    for (j in 1:length(newstates))
        for (k in newstates[[j]])
            index[k] <- j
    state.names <- rownames(object$trans)
    if (is.null(state.names)) state.names <- 1:nrow(object$trans)
    new.state.names <- tapply(state.names, index, paste0, collapse=sep)
    trans <- object$trans
    for (i in rev(which[-1])) {
        trans <- trans[-i,,drop=FALSE]
        trans <- trans[,-i,drop=FALSE]
    }
    ## sum <- function(x) tapply(x, index, base::sum)
    sum3 <- function(x,index,DIM=2) {
        dims <- dim(x)
        dims[DIM] <- length(unique(index))
        y <- array(0,dims)
        dimnames <- dimnames(x)
        dimnames[[DIM]] <- new.state.names
        for (i in 1:length(index)) {
            if (DIM==1)
                y[index[i],,] <- y[index[i],,] + x[i,,]
            else if (DIM==2)
                y[,index[i],] <- y[,index[i],] + x[,i,]
        }
        dimnames(y) <- dimnames
        y
    }
    sum4 <- function(x,index) {
        ## assumes DIM=2
        dims <- dim(x)
        dims[2] <- length(unique(index))
        dimnames <- dimnames(x)
        dimnames[[2]] <- new.state.names
        y <- array(0,dims)
        for (i in 1:length(index))
            y[,index[i],,] <- y[,index[i],,] + x[,i,,]
        dimnames(y) <- dimnames
        y
    }
    P <- sum3(object$P,index)
    L <- sum3(object$L,index)
    Pu <- sum4(object$Pu,index)
    Lu <- sum4(object$Lu,index)
    if ("C" %in% names(object)) {
        C <- sum3(object$C,index)
        Cu <- sum4(object$Cu,index)
    }
    index <- 1
    for (i in 1:nrow(object$trans)) {
        if (any(i==trans,na.rm=TRUE)) {
            trans[which(i==trans)] <- index
            index <- index + 1
        }
    }
    colnames(trans) <- rownames(trans) <- new.state.names
    newobject <- structure(list(time = object$time,
                                P=P, L=L, Pu=Pu, Lu=Lu,
                                res=NULL, vcov=object$vcov,
                                newdata=object$newdata,
                                trans=trans, call=call),
                           class=class(object))
    newobject
}

stratifiedModel <- function (object,strata) {
  ## if(inherits (object, "stratifiedModel"))
  ##     return(object)
  if(! inherits (object, "stratifiedModel"))
    class (object) <- c("stratifiedModel", class (object))
  object$strata.name <- substitute (strata)
  stopifnot(is.name(object$strata.name))
  object$strata.index <- NULL
  object
}
"[[.stratifiedModel" <- function (object, index) {
  object$strata.index <- index
  object
}
## does the following only work for index 1..K?
predict.stratifiedModel <- function (object, type, newdata) {
    if (! is.null(object$strata.index))
        newdata[[object$strata.name]] <- object$strata.index
    ## otherwise assume that the strata is specified in newdata
    NextMethod ("predict", object)
}
## predict.test <- function(object, type, newdata) return(newdata$strata)
## m <- stratifiedModel("class<-"(list(),"test"), strata)
## predict(m[[2]],"",data.frame(x=1))
coef.stratifiedModel <- function(object, ...) NextMethod("coef", object)
vcov.stratifiedModel <- function(object, ...) NextMethod("vcov", object)
