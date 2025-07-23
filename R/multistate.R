markov_msm <- 
    function (x, trans, t = c(0,1), newdata = NULL, init=NULL,
              tmvar = NULL, 
              sing.inf = 1e+10, method="adams", rtol=1e-10, atol=1e-10, slow=FALSE,
              min.tm=1e-8,
              utility=function(t) rep(1, nrow(trans)),
              utility.sd=rep(0,nrow(trans)),
              use.costs=FALSE,
              transition.costs=function(t) rep(0,sum(!is.na(trans))), # per transition
              transition.costs.sd=rep(0,sum(!is.na(trans))),
              state.costs=function(t) rep(0,nrow(trans)), # per unit time
              state.costs.sd=rep(0,nrow(trans)),
              discount.rate = 0,
              block.size=500,
              spline.interpolation=FALSE,
              debug = FALSE,
              ...)
{
    stopifnot(requireNamespace("deSolve", quietly=TRUE))
    call <- match.call()
    inherits <- function(x, ...)
        base::inherits(x, ...) ||
            (base::inherits(x, c("hazFun","zeroModel","hrModel","stratifiedModel"))
                && base::inherits(x$base, ...))
    base.classes <- c("stpm2","pstpm2","glm","survPen","gam","aft","flexsurvreg","aftreg","smoothpwc")
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
                                              transition.costs.sd = transition.costs.sd,
                                              state.costs.sd = state.costs.sd,
                                              state.costs = state.costs,
                                              discount.rate = discount.rate,
                                              block.size = block.size,
                                              spline.interpolation = spline.interpolation,
                                              ...))
        return(do.call(rbind.markov_msm,lst))
    }
    if (use.costs)
        slow <- TRUE
    nt <- length(t)
    if (nt < 2) 
        stop("number of times should be at least two")
    stopifnot(length(utility(t[2])) %in% c(1,nrow(trans)))
    if (is.null(tmvar) && all(sapply(x,inherits,c("stpm2","pstpm2","aft","survPen","flexsurvreg",
                                                  "aftreg","smoothpwc"))))
        tmvar <- sapply(x,function(object)
            if(inherits(object,c("stpm2","pstpm2"))) object@timeVar
            else if (inherits(object,"aft")) object@args$timeVar
            else if (inherits(object,"flexsurvreg")) "t"
            else if (inherits(object,"aftreg"))
                local({lhs <- object$call$formula[[2]]
                    deparse(if(length(lhs)==4) lhs[[4]] else lhs[[3]])})
            else if (inherits(object,"smoothpwc")) "t"
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
    if (spline.interpolation) {
        replace <- function(object,tmvar) {
            logModels <- c("aftreg","flexsurvreg","stpm2","pstpm2","aft")
            logModel <- inherits(object, logModels)
            t <- if (logModel) exp(seq(log(max(c(min.tm,min(t)))), log(max(t)), length.out=1000))
                 else seq(max(c(min.tm,min(t))), max(t), length.out=1000)
            if (inherits(object,"hazFun")) object
            else if (inherits(object, "addModel"))
                structure(lapply(object,replace,tmvar=tmvar), class="addModel")
            else makeSplineFun(object, newdata=newdata, tmvar=tmvar, min.tm=min.tm, tm=t,
                               log=logModel)
        }
        x <- mapply(replace, x, tmvar, SIMPLIFY=FALSE)
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
        unlist(lapply(objects,function(x) coef(x))) # sapply(objects, coef) FAILS
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
        Q <- array(0, c(nstates, nstates, nobs))
        Qu <- array(0, c(nstates, nstates, ncoef, nobs))
        QC <- array(0, c(nstates, nobs))
        QuC <- array(0, c(nstates, ncoef, nobs))
        PQCu <- array(0, c(nstates, ntr, nobs))
        froms <- tos <- rep(0,ntr)
        for (i in 1:ntr) {
            newdata[[tmvar[i]]] <- pmax(min.tm,t) # special case: t=0
            ij <- matwhich(trans == i)
            froms[i] <- from <- ij[1]
            tos[i] <- to <- ij[2]
            Q[from,to,] <- pmax(0,predict(x[[i]], newdata=newdata, type="haz", tmvar=tmvar[i]))
            Qu[from,to,1:cs[i]+cumcs[i],] <- t(predict(x[[i]], newdata=newdata, type="gradh",
                                                       tmvar=tmvar[i])) # transposed
            if (use.costs) {
                ## collapse across "to" states
                QC[from,] <- QC[from,] + Q[from,to,]*transition.costs(t)[[i]]
                QuC[from,1:cs[i]+cumcs[i],] <- QuC[from,1:cs[i]+cumcs[i],] + Qu[from,to,1:cs[i]+cumcs[i],]*transition.costs(t)[i]
                PQCu[from,i,] <- P[from,]*Q[from,to,]
            }
        }
        Q[is.na(Q)] <- Qu[is.na(Qu)] <- QC[is.na(QC)] <- QuC[is.na(QuC)] <- 0
        Q[is.infinite(Q) & Q > 0] <- Qu[is.infinite(Qu) & Qu > 0] <- sing.inf
        QC[is.infinite(QC) & QC > 0] <- QuC[is.infinite(QuC) & QuC > 0] <- sing.inf
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
                    QCi <- QC[,i]                # nstates
                    QuCi <- QuC[,,i]             # nstates*ncoef
                    dTCdt[,i] <- Pi*(state.costs(t) + QCi)      # nstates
                    dTCudt[,,i] <- Pui*(state.costs(t) + QCi) +
                        Pi*QuCi # nstates*ncoef
                }
            }
        } else {
            Qu2 <- array(Qu,c(nstates*nstates,ncoef,nobs))
            QuC2 <- array(QuC,c(nstates*nstates,ncoef,nobs))
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
                   P*exp(-lambda.discount*t), # d/dt for gradient of U wrt utilities *and* gradient of C wrt state costs
                   dTCdt*exp(-lambda.discount*t),
                   dTCudt*exp(-lambda.discount*t),
                   PQCu*exp(-lambda.discount*t) # gradient of C wrt transition costs
                   ))
        else
            list(c(dPdt,
                   dPudt,
                   P*utility(t)*exp(-lambda.discount*t),
                   Pu*utility(t)*exp(-lambda.discount*t),
                   P*exp(-lambda.discount*t) # d/dt for gradient of U wrt utilities (assumes u'_m(t)=1)
                   ))
    }
    class(x) <- "surv_list"
    Sigma <- vcov(x)
    if (any(is.na(Sigma)))
        stop("NAs in the covariance matrix")
    if (any(is.na(coef(x))))
        stop("NAs in the coefficient vector")
    init <- if (use.costs) c(rep(init,nobs),
                             rep(0,3*nobs*nstates+3*nobs*nstates*ncoef+nobs*nstates*ntr))
            else
                c(rep(init,nobs),rep(0,2*nobs*nstates+2*nobs*nstates*ncoef))
    res <- deSolve::ode(y = init, times = t, func = dp, method=method, atol=atol, rtol=rtol,
               ...)
    P <- array(res[,1+1:(nstates*nobs)],c(nt,nstates,nobs))
    Pu <- array(res[,1+nobs*nstates+1:(nobs*nstates*ncoef)],c(nt,nstates,ncoef,nobs))
    L <- array(res[,1+nstates*nobs+nstates*nobs*ncoef+1:(nstates*nobs)],c(nt,nstates,nobs))
    Lu <- array(res[,1+2*nstates*nobs+nstates*nobs*ncoef+1:(nobs*nstates*ncoef)],c(nt,nstates,ncoef,nobs))
    Pdisc <- array(res[,1+2*nstates*nobs+2*nstates*nobs*ncoef+1:(nstates*nobs)],c(nt,nstates,1,nobs))
    state.names <- if(is.null(rownames(trans))) 1:nrow(trans) else rownames(trans)
    dimnames(P) <- dimnames(L) <- 
        list(time=t, state= state.names, obs=rownames(newdata))
    dimnames(Pdisc) <- 
        list(time=t, state= state.names, coef="Pdisc", obs=rownames(newdata))
    coef.names <- unlist(lapply(x, function(xi) names(coef(xi))))
    dimnames(Pu) <- dimnames(Lu) <-
        list(time=t,
             state=state.names,
             coef=coef.names,
             obs=rownames(newdata))
    if (use.costs) {
        costs <- array(res[,1+3*nstates*nobs+2*nstates*nobs*ncoef+1:(nstates*nobs)],c(nt,nstates,nobs))
        costsu <- array(res[,1+4*nstates*nobs+2*nstates*nobs*ncoef+1:(nobs*nstates*ncoef)],c(nt,nstates,ncoef,nobs))
        costsu.trans <-   array(res[,1+4*nstates*nobs+3*nstates*nobs*ncoef+1:(nobs*nstates*ntr)],c(nt,nstates,ntr,nobs))
        dimnames(costs) <- dimnames(P)
        dimnames(costsu) <- dimnames(Pu)
        dimnames(costsu.trans) <- list(time=t,
                                       state=state.names,
                                       trans=rownames(trans),
                                       obs=rownames(newdata))
    } 
    structure(list(time = t, P=P, Pu=Pu, L=L, Lu=Lu, Pdisc=Pdisc,
                   res=if (debug) res else NULL,
                   use.costs=use.costs,
                   costs=if (use.costs) costs else NULL,
                   costsu=if (use.costs) costsu else NULL,
                   costsu.trans=if (use.costs) costsu.trans else NULL,
                   vcov=Sigma, newdata=newdata, trans=trans,
                   utility.sd=utility.sd,
                   state.costs.sd=state.costs.sd,
                   transition.costs.sd=transition.costs.sd,
                   coefficients=coef(x),
                   call=call, x=x),
              class="markov_msm")
}

## Pu, vcov
## Lu || Pdisc, bdiag(vcov,diag(utility.sd^2))
## costsu || Pdisc || costsu.trans, bdiag(vcov, diag(state.costs.sd^2), diag(transition.costs.sd^2))

bdiag <- function(...,initial=0) {
    l <- list(...)
    for (elt in l)
        if (is.array(elt) && length(dim(elt))>2)
            stop("not defined for arrays with more than 2 dimensions")
    l <- lapply(l, as.matrix)
    nr <- sapply(l,nrow)
    nc <- sapply(l,ncol)
    cumnr <- c(0,cumsum(nr))
    cumnc <- c(0,cumsum(nc))
    out <- matrix(initial,sum(nr),sum(nc))
    for (i in 1:length(l))
        out[(cumnr[i]+1):cumnr[i+1], (cumnc[i]+1):cumnc[i+1]] <- l[[i]]
    ## if rownames/colnames are all missing, then NULL else convert the NULLs to blank strings
    replace <- function(l,f,len) {
        lf <- lapply(l,f)
        nulls <- sapply(lf, is.null)
        if (all(nulls)) NULL
        else unlist(lapply(1:length(l), function(i) if (nulls[i]) rep("",len(l[[i]])) else lf[[i]]))
    }
    rownames(out) <- replace(l,rownames, nrow)
    colnames(out) <- replace(l,colnames, ncol)
    out
}
## all(rownames(bdiag(matrix(1:4,2),c(a=5,b=6),7L,t(8:9),initial=0L))==c("","","a","b","",""))
## is.integer(bdiag(matrix(1:4,2),a=5:6,7L,t(8:9),initial=0L))
## bdiag(list(1:2),list(3:4+0),initial=list())

hazFun <- function(f, tmvar="t", ...) {
    nms <- formalArgs(f)
    if(all(c("newdata","t") %in% nms))
        newf <- function(newdata) f(t=newdata[[tmvar]],newdata=newdata,...)
    else if("t" %in% nms)
        newf <- function(newdata) f(t=newdata[[tmvar]], ...)
    else if("newdata" %in% nms)
        newf <- function(newdata) f(newdata=newdata, ...)
    else 
        newf <- function(newdata) f(...)
    out <- list(haz=newf)
    class(out) <- "hazFun"
    out
}
predict.hazFun <- function(object, newdata, type=c("haz","gradh"), ...) {
    type <- match.arg(type)
    val <- if (type=="haz") object$haz(newdata) else 0
    if (length(val)==1 && length(val)<nrow(newdata))
        val <- rep(val,nrow(newdata))
    val
}
coef.hazFun <- function(object, ...) c(hazFun=0)
vcov.hazFun <- function(object, ...) matrix(0,1,1,FALSE,list("hazFun","hazFun"))


splinefunx <- function(x, y, method="natural", constant.left=FALSE, constant.right=FALSE, ...) {
    xstar <- x
    ystar <- y
    if (constant.right) {
        xstar <- c(xstar,tail(x,1)+(1:10)*1e-7)
        ystar <- c(ystar,rep(tail(y,1),10))
    }
    if (constant.left) {
        xstar <- c(x[1]-(10:1)*1e-7,xstar)
        ystar <- c(rep(y[1],10),ystar)
    }
    stats::splinefun(xstar, ystar, ..., method=method)
}
smoothpwc <- function(midts, rates, tmvar="t", offsetvar="", ...) {
    log.smoother <- splinefunx(midts, log(rates), constant.right=TRUE)
    haz <- function(newdata) {
        t <- newdata[[tmvar]] + (if (offsetvar!="") newdata[[offsetvar]] else 0)
        exp(log.smoother(t))
    }
    structure(list(haz=haz), class="smoothpwc")
}
coef.smoothpwc <- function(object, ...) c(smoothpwc=0)
vcov.smoothpwc <- function(object, ...) matrix(0,1,1,FALSE,list("pwc","pwc"))
predict.smoothpwc <- function(object, newdata, type=c("haz","gradh"), ...) {
    type <- match.arg(type)
    val <- if (type=="haz") object$haz(newdata) else 0
    if (length(val)==1 && length(val)<nrow(newdata))
        val <- rep(val,nrow(newdata))
    val
}


as.data.frame.markov_msm <- function(x, row.names=NULL, optional=FALSE,
                                     ci=TRUE,
                                     P.conf.type="logit", L.conf.type="log", C.conf.type="log",
                                     P.range=c(0,1), L.range=c(0,Inf),
                                     C.range=c(0,Inf),
                                     state.weights=NULL,
                                     obs.weights=NULL,
                                     ...) {
    if(!is.null(obs.weights))
        x <- standardise(x, weights=obs.weights, ...)
    if (!is.null(state.weights))
        return(state_weights(x=x, row.names=row.names,
                             optional=optional, ci=ci,
                             P.conf.type=P.conf.type, L.conf.type=L.conf.type,
                             C.conf.type=C.conf.type,
                             P.range=P.range, L.range=L.range, C.range=C.range,
                             state.weights=state.weights, ...))
    state.names <- rownames(x$trans)
    perm <- function(x) aperm(x, c(2,3,1))
    ## probabilities
    seP <- apply(x$Pu, c(1,2), function(m) sqrt(colSums(m * (vcov(x) %*% m)))) # c(nobs,nt,nstates)
    seP <- array(seP, dim(x$Pu)[c(4,1,2)]) # c(nobs,nt,nstates))
    seP <- aperm(seP,c(2,3,1))
    ## utilities
    keep.ith <- function(a,i) {
        slice <- a[,i,,] 
        a <- a*0
        a[,i,,] <- slice
        a
    }
    Lu <- do.call(abind, c(list(3, x$Lu),
                           lapply(1:nrow(x$trans), function(i) keep.ith(x$Pdisc,i))))
    vcovL <- bdiag(vcov(x), diag(x$utility.sd^2))
    seL <- apply(Lu, c(1,2), function(m) sqrt(colSums(m * (vcovL %*% m)))) # c(nobs,nt,nstates)
    seL <- array(seL, dim(Lu)[c(4,1,2)]) # c(nobs,nt,nstates))
    seL <- aperm(seL,c(2,3,1))
    out <- data.frame(P=as.vector(perm(x$P)),
                      seP=as.vector(perm(seP)),
                      L=as.vector(perm(x$L)),
                      seL=as.vector(perm(seL)))
    if (!is.null(x$costsu)) {
        ## error in the following 
        costsu <- do.call(abind, c(list(3,x$costsu),
                                   lapply(1:nrow(x$trans), function(i) keep.ith(x$Pdisc,i)),
                                   list(x$costsu.trans)))
        vcovC <- bdiag(vcov(x), diag(x$state.costs.sd^2), diag(x$transition.costs.sd^2))
        seC <- apply(costsu, c(1,2), function(m) sqrt(colSums(m * (vcovC %*% m)))) # c(nobs,nt,nstates)
        seC <- array(seC, dim(costsu)[c(4,1,2)]) # c(nobs,nt,nstates))
        seC <- aperm(seC,c(2,3,1))
        out$C <- as.vector(perm(x$costs))
        out$seC <- as.vector(perm(seC))
    }
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
        if (!is.null(x$costs)) {
            tmp <- surv.confint(out$C,out$seC, conf.type=C.conf.type, min.value=C.range[1], max.value=C.range[2])
            out$C.lower=tmp$lower
            out$C.upper=tmp$upper
        }
    }
    newdata <- as.data.frame(x$newdata[as.numeric(out$obs), , drop=FALSE])
    names(newdata) <- names(x$newdata)
    out <- "rownames<-"(cbind(newdata, out),rownames(out))
    if(!is.null(row.names)) rownames(out) <- row.names
    out
}

# if var(P1)=g1*Sigma*g1' and var(P2)=g2*Sigma*g2' then var(P1+P2)=(g1+g2)*Sigma*(g1+g2)' 
state_weights <- function(x, row.names=NULL, optional=FALSE,
                          ci=TRUE,
                          P.conf.type="logit", L.conf.type="log", C.conf.type="log",
                          P.range=c(0,1), L.range=c(0,Inf),
                          C.range=c(0,Inf),
                          state.weights,
                          ...) {
    stopifnot(!missing(state.weights))
    state.names <- rownames(x$trans)
    ## stopifnot(length(state.names) == length(state.weights))
    w <- function(m) colSums(m * state.weights)
    perm <- function(x) aperm(x, c(2,3,1))
    ## probabilities
    seP <- apply(x$Pu, c(1,4), function(m) sqrt(colSums(w(m) * (vcov(x) %*% w(m))))) # c(nobs,nt)
    ## utilities
    keep.ith <- function(a,i) {
        slice <- a[,i,,] 
        a <- a*0
        a[,i,,] <- slice
        a
    }
    Lu <- do.call(abind, c(list(3, x$Lu),
                           lapply(1:nrow(x$trans), function(i) keep.ith(x$Pdisc,i))))
    vcovL <- bdiag(vcov(x), diag(x$utility.sd^2))
    seL <- apply(Lu, c(1,4), function(m) sqrt(colSums(w(m) * (vcovL %*% w(m))))) # c(nobs,nt)
    out <- data.frame(P=as.vector(apply(x$P,1,w)),
                      seP=as.vector(t(seP)),
                      L=as.vector(apply(x$L,1,w)),
                      seL=as.vector(t(seL)))
    if (!is.null(x$costsu)) {
        ## error in the following 
        costsu <- do.call(abind, c(list(3,x$costsu),
                                   lapply(1:nrow(x$trans), function(i) keep.ith(x$Pdisc,i)),
                                   list(x$costsu.trans)))
        vcovC <- bdiag(vcov(x), diag(x$state.costs.sd^2), diag(x$transition.costs.sd^2))
        seC <- apply(costsu, c(1,4), function(m) sqrt(colSums(w(m) * (vcovC %*% w(m))))) # c(nobs,nt)
        out$C <- as.vector(apply(x$costs,1,w))
        out$seC <- as.vector(t(seC))
    }
    dimnames. <- dimnames(x$P)[c(3,1)] 
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
        if (!is.null(x$costs)) {
            tmp <- surv.confint(out$C,out$seC, conf.type=C.conf.type, min.value=C.range[1], max.value=C.range[2])
            out$C.lower=tmp$lower
            out$C.upper=tmp$upper
        }
    }
    newdata <- as.data.frame(x$newdata[as.numeric(out$obs), , drop=FALSE])
    names(newdata) <- names(x$newdata)
    out <- "rownames<-"(cbind(newdata, out),rownames(out))
    if(!is.null(row.names)) rownames(out) <- row.names
    out
}


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
    Pdisc <- add.dim(apply(x$Pdisc,1:3,wtd.mean))
    dimnames(P) <- dimnames(L) <- list(time=x$time, state=state.names, obs=1)
    oldnames <- dimnames(x$Pu)
    dimnames(Pdisc) <- list(time=x$time, state=state.names, coef="Pdisc", obs=1)
    dimnames(Pu) <- dimnames(Lu) <-
        list(time=x$time, state=state.names, coef=oldnames[[3]], obs=1)
    ## update x
    x$newdata <- x$newdata[0,,drop=FALSE]
    x$P <- P
    x$Pu <- Pu
    x$L <- L
    x$Lu <- Lu
    x$Pdisc <- Pdisc
    x$call <- call
    x
}
vcov.survPen <- function(object) object$Vp
"coef<-.survPen" <- function(x,value) {
    x$coefficients <- value
    x
}
survPenWrap <- function(object) {
    stopifnot(inherits(object, "survPen"))
    if (!inherits(object, "survPenWrap"))
        class(object) <- c("survPenWrap", class(object))
    object
}
predict.survPenWrap <- function(object, newdata, type=NULL, min.eps=1e-6, ...) {
    newdata[[object$t1.name]] <- pmax(min.eps, newdata[[object$t1.name]])
    if (is.null(type)) NextMethod("predict", object)
    else if (type=="haz") predict(object, newdata=newdata, do.surv=FALSE)$haz
    else if (type=="gradh") predict(object, newdata=newdata, type="lpmatrix") *
                                predict(object, newdata=newdata, do.surv=FALSE)$haz
    else NextMethod("predict", object)
}
predict.glm <- function (object, newdata=NULL, type=NULL, ...) {
    if (is.null(type)) NextMethod("predict", object)
    if(type=="gradh" && !(object$family$link %in% c("identity","log")))
        stop("Currently only implemented for log and identity links")
    ## stopifnot() # Poisson family with log link?
    ## assumes response is a rate
    if(!is.na(pmatch(type,"hazard")))
        stats::predict.glm(object, newdata=newdata, type="response", ...)
    else if (type=="gradh" && object$family$link == "log")
        stats::predict.glm(object, newdata=newdata, type="response", ...) *
            lpmatrix.lm(object, newdata=newdata)
    else if (type=="gradh" && object$family$link == "identity")
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
    if (type=="haz") predict(object, newdata=newdata, type="response", ...) 
    else if (type=="gradh") {
        as.vector(predict(object, newdata=newdata, type="response", ...)) *
            predict(object, newdata=newdata, type="lpmatrix", ...)
    } else NextMethod("predict", object)
}
predict.flexsurvreg <- function(object, newdata, type=NULL, t=NULL, se.fit=FALSE, tmvar="t", ...) {
    if (is.null(t)) t <- unique(newdata[[tmvar]])
    if (is.null(type)) return(summary(object, newdata=newdata, t=t, ci=se.fit, ...))
    if (type=="haz")
        sapply(summary(object, newdata=newdata, t=t, ci=se.fit, type="haz"),
               "[[", "est")
    else if (type=="gradh") {
        g <- function(coef) {
            object$res[,"est"] <- object$res.t[,"est"] <- coef
            predict(object, newdata, type="haz", t=t, ci=FALSE, ...)
        }
        ## numDeriv::jacobian(g, x=coef(object), method="simple")
        grad(g, coef(object))
    }
    else summary(object, newdata=newdata, type=type, t=t, ci=se.fit, ...)
}

vcov.aftreg <- function(object, ...)
    object$var

predict.aftreg <- function (object, type = c("haz", "cumhaz", "density", "surv", "gradh"), 
                            newdata, t=NULL, tmvar=NULL, na.action=na.pass, fd=FALSE, ...) 
{
    ## utility function
    mref <- function(m,i,j) i+(j-1)*nrow(m)
    if (!inherits(object, "aftreg")) 
        stop("Works only with 'aftreg' objects.")
    type <- match.arg(type)
    if (type=="gradh" && fd) {
        coef <- coef(object)
        g1 <- function(coef) {
            object$coefficients <- coef
            predict.aftreg(object, newdata=newdata, type="haz", t=t,
                           tmvar=tmvar, na.action=na.action, fd=FALSE, ...)
        }
        ## return(numDeriv::jacobian(g, coef))
        ## return(numDeriv::jacobian(g, coef, method="simple"))
        return(grad(g1, coef))
    }
    if (is.null(t)) {
        if (is.null(tmvar)) {
            lhs <- object$call$formula[[2]]
            expr <- if (length(lhs)==4) lhs[[3]] else lhs[[2]]
            t <- eval(expr, newdata, parent.frame())
            rm(lhs,expr)
        } else
            t <- newdata[[tmvar]]
    }
    if (length(t)<nrow(newdata)) t <- rep(t,length=nrow(newdata))
    if (type=="haz" && fd) {
        S <- predict.aftreg(object, newdata=newdata, type="surv", t=t,
                            tmvar=tmvar, na.action=na.action, ...)
        g2 <- function(t)
            predict.aftreg(object, newdata=newdata, type="surv", t=t,
                           tmvar=tmvar, na.action=na.action, fd=FALSE, ...)
        ## return(-numDeriv::jacobian(g, t)/S)
        return(-grad(g2, t)/S)
    }
    ## BUG: add names to object$levels
    if (any(object$isF))
        names(object$levels) <- object$covars[object$isF]
    Terms <- object$terms
    if (!inherits(Terms, "terms")) 
        stop("invalid terms component of  object")
    strats <- attr(Terms, "specials")$strata
    intercept <- attr(Terms, "intercept")
    Terms <- delete.response(Terms)
    newframe <- stats::model.frame(Terms, data = newdata, 
                                   na.action = na.action,
                                   xlev = object$levels)
    if (length(strats)) {
        if (!("strata" %in% names(object)))
            stop("predictions with strata requires 'aftreg(..., x=TRUE)'")
        temp <- untangle.specials(Terms, "strata", 1)
        if (length(temp$vars) == 1) 
            strata.keep <- newframe[[temp$vars]]
        else strata.keep <- survival::strata(newframe[, temp$vars], shortlabel = TRUE)
        strata.keep <- factor(as.character(strata.keep),levels=levels(object$strata))
        strata <- as.numeric(strata.keep)
        ## remove strats from object$terms
        newTerms <- Terms[-temp$terms]
        attr(newTerms, "intercept") <- attr(Terms, "intercept")
        Terms <- newTerms
    } else {
        stopifnot(object$n.stata == 1)
        strata <- rep(1, nrow(newframe))
    }
    "%mv%" <- function(a,b) as.vector(a %*% b)
    newframe <- stats::model.frame(Terms, data = newdata, 
                                   na.action = na.action,
                                   xlev = object$levels)
    object$terms <- Terms
    ncov <- length(object$means)
    if (object$pfixed) {
        p <- object$shape[strata]
        lambda <- exp(object$coefficients[ncov + strata])
    } else  {
        p <- exp(object$coefficients[ncov + strata * 2])
        lambda <- exp(object$coefficients[ncov + strata * 2 - 1])
    }
    if (ncov) {
        x <- model.matrix(object, newframe)[,-1,drop=FALSE] # is the intercept always first?
        ## x <- sweep(x, 2, object$means) # no longer sweep:)
        param.scale <- if (object$param=="lifeAcc") -1 else 1
        lambda <- lambda *
            exp(param.scale*(x %mv% object$coefficients[1:ncov]))
    } else {
        x <- NULL
    }
    xx <- t
    if (object$dist == "weibull") {
        dist <- "Weibull"
        haza <- eha::hweibull
        Haza <- eha::Hweibull
        Surviv <- stats::pweibull
        Dens <- stats::dweibull
    }
    else if (object$dist == "loglogistic") {
        dist <- "Loglogistic"
        haza <- eha::hllogis
        Haza <- eha::Hllogis
        Surviv <- eha::pllogis
        Dens <- eha::dllogis
    }
    else if (object$dist == "lognormal") {
        dist = "Lognormal"
        haza <- eha::hlnorm
        Haza <- eha::Hlnorm
        Surviv <- stats::plnorm
        Dens <- stats::dlnorm
    }
    else if (object$dist == "ev") {
        dist = "Extreme value"
        haza <- eha::hEV
        Haza <- eha::HEV
        Surviv <- eha::pEV
        Dens <- eha::dEV
    }
    else if (object$dist == "gompertz") {
        dist = "Gompertz"
        haza <- eha::hgompertz
        Haza <- eha::Hgompertz
        Surviv <- eha::pgompertz
        Dens <- eha::dgompertz
        ## canonical parameterisation
        p <- p/lambda
    }
    if (type=="haz")
        return(haza(xx, scale = lambda, shape = p))
    if (type=="cumhaz")
        return(Haza(xx, scale = lambda, shape = p))
    if (type=="density") {
        if (object$dist == "lognormal") {
            sdlog <- 1/p
            meanlog <- log(lambda)
            den <- stats::dlnorm(xx, meanlog, sdlog)
        }
        else 
            den <- Dens(xx, scale = lambda, shape = p)
        return(den)
    }
    if (type=="surv") {
        if (object$dist == "lognormal") {
            sdlog <- 1/p
            meanlog <- log(lambda)
            sur <- stats::plnorm(xx, meanlog, sdlog, lower.tail = FALSE)
        }
        else {
            sur <- Surviv(xx, scale = lambda, shape = p, 
                          lower.tail = FALSE)
        }
        return(sur)
    }
    if (type=="gradh" && !fd) {
        mapper <- list(ev="extreme",
                       weibull="weibull",
                       loglogistic="loglogistic",
                       Lognormal="lognormal")
        if (object$dist %in% names(mapper)) {
            distname <- mapper[[object$dist]]
            dd <- survival::survreg.distributions[[distname]]
            if (is.null(dd$itrans)) {
                trans <- function(x) x
                itrans <- function(x) x
                dtrans <- function(x) 1
            }
            else {
                trans <- dd$trans
                itrans <- dd$itrans
                dtrans <- dd$dtrans
            }
            if (!is.null(dd$dist)) 
                dd <- survival::survreg.distributions[[dd$dist]]
            pred <- log(lambda)
            scale <- 1/p
            u <- (trans(t)-pred) / scale # check dimensions
            density <- dd$density(u, list()) # F, 1-F, f, f'/f, f''/f
            S <- density[,2]
            f <- density[,3]
            f.prime <- density[,4]*f
            ## NOTE: x may be NULL
            grad.beta <- cbind(x,-1)*(f.prime*S+f^2)*dtrans(t)/(S*scale)^2
            out <- matrix(0,nrow(newdata),length(coef(object)))
            nc <- ncol(grad.beta)
            if (ncov)
                out[,1:(nc-1)] <- grad.beta[,1:(nc-1)]
            out[mref(out,1:nrow(newdata),nc+strata*2-2)] <- grad.beta[,nc]
            if (!object$pfixed) {
                grad.logscale <- (-f.prime*u*S/dtrans(t) -
                                  f^2*u/dtrans(t) -
                                    f*S/dtrans(t)) / (S*scale/dtrans(t))^2*scale
                out[mref(out,1:nrow(newdata),nc+strata*2-1)] <- -grad.logscale
            }
            return(out)
        } else {
            ## gompertz
            h <- haza(xx, scale = lambda, shape = p)
            grad.beta <- cbind(x,-1)*h*xx/p
            grad.logscale <- h
            out <- matrix(0,nrow(newdata),length(coef(object)))
            nc <- ncol(grad.beta)
            if (!is.null(x))
                out[,1:(nc-1)] <- grad.beta*x
            out[mref(out,1:nrow(newdata),nc+strata*2-2)] <- grad.beta[,nc]
            out[mref(out,1:nrow(newdata),nc+strata*2-1)] <- -grad.logscale
            return(out)
        }
    }
}

surv.confint <- function(p, se, conf.type=c("log","log-log","plain","logit","arcsin"),
                         conf.int=0.95, min.value=0, max.value=1) {
    stopifnot(length(p) == length(se))
    conf.type <- match.arg(conf.type)
    zval <- qnorm(1 - (1 - conf.int)/2)
    selog <- ifelse(p==0, NaN, se/p) # transform to selog
    out <- if (conf.type == "plain") {
               se2 <- se * zval
               list(lower = pmax(p - se2, min.value), upper = pmin(p + se2, max.value))
           }
           else if (conf.type == "log") {
               xx <- ifelse(p == 0, NaN, p)
               selog2 <- zval * selog
               temp1 <- exp(log(xx) - selog2)
               temp2 <- exp(log(xx) + selog2)
               list(lower = temp1, upper = pmin(temp2, max.value))
           }
           else if (conf.type == "log-log") {
               xx <- ifelse(p == 0, NaN, p)
               selog2 <- zval * selog/log(xx)
               temp1 <- exp(-exp(log(-log(xx)) - selog2))
               temp2 <- exp(-exp(log(-log(xx)) + selog2))
               list(lower = temp1, upper = temp2)
           }
           else if (conf.type == "logit") {
               xx <- ifelse(p == 0 | p == 1, NaN, p)
               selog2 <- zval * selog * (1 + xx/(1 - xx))
               temp1 <- 1 - 1/(1 + exp(log(xx/(1 - xx)) - selog2))
               temp2 <- 1 - 1/(1 + exp(log(xx/(1 - xx)) + selog2))
               list(lower = temp1, upper = temp2)
           }
           else if (conf.type == "arcsin") {
               xx <- ifelse(p == 1, NaN, p)
               selog2 <- 0.5 * zval * selog * sqrt(xx/(1 - xx))
               list(lower = (sin(pmax(0, asin(sqrt(xx)) - selog2)))^2,
                    upper = (sin(pmin(pi/2, asin(sqrt(xx)) + selog2)))^2)
           }
    if (any(index <- is.na(p) | is.na(se)))
        out$lower[index] <- out$upper[index] <- NA
    out
}
print.markov_msm <- function(x, 
                             digits=5,
                             se=FALSE, ci=TRUE,
                             P.conf.type="logit", L.conf.type="log",
                             C.conf.type="log",
                             P.range=c(0,1), L.range=c(0,Inf),
                             C.range=c(0,Inf),
                             ...) {
    df <- as.data.frame(x, ci=ci,
                        P.conf.type=P.conf.type, L.conf.type=L.conf.type,
                        C.conf.type=C.conf.type,
                        P.range=P.range, L.range=L.range, C.range=C.range)
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
    0.0*predict(object$base, type=type, ...) 
}
## we keep the same length/dim to support comparisons between different interventions
coef.zeroModel <- function(object, ...) 0.0*coef(object$base, ...)
vcov.zeroModel <- function(object, ...) 0.0*vcov(object$base, ...)

hrModel <- function(object, hr=1, ci=NULL, seloghr=NULL) {
    stopifnot(is.null(ci) || (is.numeric(ci) && length(ci)==2))
    stopifnot(is.numeric(hr) && length(hr)==1 && hr>0)
    newobject <- if (inherits(object, "hrModel")) object else structure(list(base=object),
                                                                        class="hrModel")
    if (!is.null(ci) && is.null(seloghr))
        seloghr <- log(ci[2]/ci[1])/2/qnorm(0.975)
    attr(newobject, "loghr") <- log(hr)
    attr(newobject, "seloghr") <- if (is.null(seloghr)) 0 else seloghr
    newobject
}
predict.hrModel <- function(object, type=c("haz","gradh"), ...) {
    type <- match.arg(type)
    hr <- exp(attr(object,"loghr"))
    pred1 <- predict(object$base, type="haz", ...)
    if (type=="haz") hr*pred1
    else
        cbind(predict(object$base, type="gradh", ...)*hr, pred1*hr)
}
## This has different lengths/dimensions to the base model:
##   wrap base intervention in hrModel(..., hr=1)
coef.hrModel <- function(object, ...) c(coef(object$base), attr(object,"loghr"))
vcov.hrModel <- function(object, ...) {
    out <- rbind(cbind(vcov(object$base),0),0)
    out[nrow(out),ncol(out)] <- attr(object, "seloghr")^2
    out
}

## aftModel requires dh/dt - performed using finite differences
## Slow? Cool for completeness - but unlikely to be used
## axiom: h := operator 'h
## axiom: hstar := h(t*exp(eta),beta)*exp(eta)
## axiom: D(hstar,beta)
## axiom: D(hstar,eta)
## h(t)=h(exp(u)) where t=exp(u)
## => dh/dt = dh/du * du/dt where du/dt=1/t
## and dh/du ~= (h(t*exp(eps)) - h(t/exp(eps)))/2/eps
aftModel <- function(object, af=1, ci=NULL, selogaf=NULL) {
    stopifnot(is.null(ci) || (is.numeric(ci) && length(ci)==2))
    stopifnot(is.numeric(af) && length(af)==1 && af>0)
    if (!is.null(ci) && is.null(selogaf))
        selogaf <- log(ci[2]/ci[1])/2/qnorm(0.975)
    newobject <- if (inherits(object, "aftModel")) object else structure(list(base=object),
                                                                        class="aftModel")
    attr(newobject, "logaf") <- log(af)
    attr(newobject, "selogaf") <- if (is.null(selogaf)) 0 else selogaf
    newobject
}
predict.aftModel <- function(object, newdata, type=c("haz","gradh"), tmvar=NULL, eps=1e-5, ...) {
    type <- match.arg(type)
    af <- exp(attr(object,"logaf"))
    time <- newdata[[tmvar]]
    timeStar <- newdata[[tmvar]] <- af*time
    pred1 <- predict(object$base, type="haz", newdata=newdata, ...)
    if (type=="haz") af*pred1
    else {
        newdata2 <- newdata
        newdata2[[tmvar]] <- timeStar*exp(eps)
        predu <- predict(object$base, type="haz", newdata=newdata2, ...)
        newdata2[[tmvar]] <- timeStar/exp(eps)
        predl <- predict(object$base, type="haz", newdata=newdata2, ...)
        dhdt <- (predu-predl)/2/eps/timeStar
        cbind(predict(object$base, type="gradh", newdata=newdata, ...)*af,
              timeStar*af*dhdt + af*pred1)
    }
}
## This has different lengths/dimensions to the base model:
##   wrap base intervention in aftModel(..., af=1)
coef.aftModel <- function(object, ...) c(coef(object$base), attr(object,"logaf"))
vcov.aftModel <- function(object, ...) {
    out <- rbind(cbind(vcov(object$base),0),0)
    out[nrow(out),ncol(out)] <- attr(object, "selogaf")^2
    out
}
## h'(t):
## gsm: S(t|x)=G(eta(t,x)) => H(t|x)=-log(S(t|x))
## => h(t|x)= -G'(eta(t,x))/G(eta(t,x))*eta_{,1}(t,x)
## => h'(t|x) = mess

## export
splineFun <- function(time, rate, method="natural", scale=1, ...) {
    fun <- stats::splinefun(time, log(rate), method=method, ...)
    hazFun(function(t) scale*exp(fun(t)))
}

## export
addModel <- function (...) {
    structure (as.list (...), class="addModel")
}
predict.addModel <- function (object, ...) {
    Reduce("+", lapply(object, predict, ...))
}
coef.addModel <- function (object)
    do.call(c,lapply(object, coef))
vcov.addModel <- function (object)
    do.call(bdiag,lapply(object, vcov))

subset.markov_msm <- function(x, subset, ...) {
    e <- substitute(subset)
    r <- eval(e, x$newdata, parent.frame())
    if (!is.logical(r)) 
        stop("'subset' must be logical")
    r <- r & !is.na(r)
    # update x
    x$P <- x$P[,,r,drop=FALSE]
    x$Pu <- x$Pu[,,,r,drop=FALSE]
    x$L <- x$L[,,r,drop=FALSE]
    x$Lu <- x$Lu[,,,r,drop=FALSE]
    x$Pdisc <- x$Pdisc[,,,r,drop=FALSE]
    x$newdata <- x$newdata[r,,drop=FALSE]
    if (x$use.costs) {
        x$costs <- x$costs[,,r,drop=FALSE]
        x$costsu <- x$costsu[,,,r,drop=FALSE]
    }
    x$res <- NULL # not strictly needed
    x$call <- match.call()
    x
}
diff.markov_msm <- function(x, y, ...) {
    stopifnot(inherits(x,"markov_msm"))
    stopifnot(inherits(y,"markov_msm"))
    ## vcov, time, trans should all be the same
    stopifnot(all(x$vcov == y$vcov))
    stopifnot(all(x$time == y$time))
    stopifnot(all(x$trans == y$trans, na.rm=TRUE))
    z <- x # copy of x
    z$P <- x$P-y$P
    z$Pu <- x$Pu-y$Pu
    z$L <- x$L-y$L
    z$Lu <- x$Lu-y$Lu
    z$Pdisc <- x$Pdisc-y$Pdisc
    z$call <- match.call()
    z$newdata[x$newdata != y$newdata] <-  NA
    if (z$use.costs) {
        z$costs <- x$costs-y$costs
        z$costsu <- x$costsu-y$costsu
    }
    class(z) <- c("markov_msm_diff","markov_msm") # not strictly "markov_msm"...
    z
}
as.data.frame.markov_msm_diff <- function(x, row.names=NULL, optional=FALSE,
                                          P.conf.type="plain", L.conf.type="plain",
                             C.conf.type="plain",
                             P.range=c(-Inf, Inf), L.range=c(-Inf, Inf),
                             C.range=c(-Inf,Inf), ...)
    as.data.frame.markov_msm(x,
                             P.conf.type=P.conf.type, L.conf.type=L.conf.type,
                             C.conf.type=C.conf.type,
                             P.range=P.range, L.range=L.range,
                             C.range=C.range, ...)

plot.markov_msm_diff <- function(x, y, ...)
    stop("Plots for markov_msm_diff have not yet been implemented")

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
    ## if (nrow(x$newdata)>1) {
    ##     warning("More than one set of covariates; covariates have been standardised")
    ##     x <- standardise(x)
    ## }
    if (ggplot2)
        ggplot.markov_msm(x, which=which, stacked=stacked, xlab=xlab, ylab=ylab,
                          alpha=alpha, ...)
    else if (lattice)
        xyplot.markov_msm(x, which=which, stacked=stacked, xlab=xlab, ylab=ylab,
                          col=col, border=border, strata=strata, ...)
    else {
        if (!missing(y)) warning("y argument is ignored")
        df <- as.data.frame(x)
        states <- unique(df$state)
        if (stacked) {
            if (which == "L")
                stop("Stacked plot is not sensible for length of stay")
            out <- graphics::plot(range(x$time),0:1, type="n", xlab=xlab, ylab=ylab, ...)
            lower <- 0
            for (i in length(states):1) { # put the last state at the bottom
                df2 <- df[df$state==states[i],]
                if (length(lower)==1) lower <- rep(0,nrow(df2))
                upper <- lower+df2[[which]]
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
                              xlab=NULL, ylab=NULL, flipped = FALSE,
                              ..., environment=parent.frame()) {
    if (requireNamespace("ggplot2", quietly=TRUE)) {
        which <- match.arg(which)
        df <- as.data.frame(data, ci=!stacked)
        if (flipped)
            df$state <- factor(df$state,levels=rev(levels(df$state)))
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
    z <- x # copy of x
    z$P <- log(x$P/y$P)
    z$L <- log(x$L/y$L) # logs of the ratios
    z$Pu <- apply(x$Pu, 3, function(slice) slice/x$P) - apply(y$Pu, 3, function(slice) slice/y$P)
    dim(z$Pu) <- dim(x$Pu)
    dimnames(z$Pu) <- dimnames(x$Pu)
    z$Lu <- apply(x$Lu, 3, function(slice) slice/x$L) - apply(y$Lu, 3, function(slice) slice/y$L)
    dim(z$Lu) <- dim(x$Lu)
    dimnames(z$Lu) <- dimnames(x$Lu)
    z$Pdisc <- log(x$Pdisc/y$Pdisc)
    z$call <- match.call()
    z$newdata[x$newdata != y$newdata] <-  NA
    class(z) <- c("markov_msm_ratio","markov_msm") # not strictly "markov_msm"...
    z
}
as.data.frame.markov_msm_ratio <- function(x, row.names=NULL, optional=FALSE, ...) {
    ## data are on a log scale!!
    z <- as.data.frame.markov_msm_diff(x, ...)
    names <- c("P","L","P.lower","P.upper","L.lower","L.upper")
    if (!is.null(x$costs))
        names <- c(names, c("C","C.lower","C.upper"))
    for (name in names)
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

abind <- function(index, ...) { # bind on nth slice for a bag of arrays
    x <- list(...)
    stopifnot(all(sapply(x,is.array)))
    dim1 <- dim(x[[1]])
    ## check index
    stopifnot(index==floor(index) && index>=1 || index<=length(dim))
    ## check all arrays have the same number of dimensions
    stopifnot(all(sapply(x[-1],function(xi) length(dim(xi))==length(dim1))))
    ## case: arrays of dimension 1
    if (length(dim1)==1) return(array(unlist(x)))
    ## check all arrays for the same dimensions (excluding the index dimension)
    for (i in ((1:length(dim1))[-index]))
        stopifnot(all(sapply(x[-1],function(xi) dim(xi)[i]==dim1[i])))
    ## re-order to put the index dimension last
    ord <- c((1:length(dim1))[-index],index)
    x <- lapply(x, function(xi) aperm(xi,ord)) # REPLACEMENT
    y <- do.call("c",lapply(x,as.vector)) # magic...
    dims <- dim1[ord]
    dimnames. <- dimnames(x[[1]])
    dims[length(dims)] <- sum(sapply(x,function(xi) dim(xi)[length(dims)]))
    if (!is.null(dimnames.[[length(dims)]]))
        dimnames.[[length(dims)]] <- unlist(lapply(x,function(xi) dimnames(xi)[length(dims)]))
    dim(y) <- dims
    dimnames(y) <- dimnames.
    aperm(y,order(ord))
}
## abind(1,array(1),array(2),array(3:4))
## abind(1,array(1,c(1,1)),array(2,c(1,1)),array(3:4,c(2,1)))
## abind(2,array(1,c(1,1)),array(2,c(1,1)),array(3:4,c(1,2)))
## abind(1,array(1,c(1,1,1)),array(2,c(1,1,1)),array(3:4,c(2,1,1)))
add_dim <- function(a,index,dimname="") {
    stopifnot(is.array(a))
    dims <- dim(a)
    ndim <- length(dims)
    dimnames. <- dimnames(a)
    stopifnot(index==floor(index) && index>=1 && index<=ndim+1)
    newdim <- if (index==1) c(1,dims)
              else if (index==ndim+1) c(dims,1)
              else c(dims[1:(index-1)], 1, dims[index:ndim])
    dimnames. <- if (index==1) c(dimname,dimnames.)
                 else if (index==ndim+1) c(dimnames.,dimname)
                 else c(dimnames.[1:(index-1)], dimname, dimnames.[index:ndim])
    array(a,newdim,dimnames.)
}
## add_dim(array(c(1,2,3,4),dimnames=list(11:14)),2,"blank")
## add_dim(array(1:4),2)
## add_dim(array(1:4,c(2,2)),1)[1,,]
## add_dim(array(1:4,c(2,2)),2)[,1,]
## add_dim(array(1:4,c(2,2)),3)[,,1]
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
    z <- x[[1]]
    z$P <- bind("P")
    z$Pu <- bind("Pu")
    z$L <- bind("L")
    z$Lu <- bind("Lu")
    z$Pdisc <- bind("Pdisc")
    z$call <- match.call()
    z$newdata <- do.call(rbind,lapply(x,"[[", "newdata"))
    ## is this a good idea? It will be lost if done more than once...
    z$newdata$.index <- unlist(lapply(1:length(x),
                                         function(i) rep(i,nrow(x[[i]]$newdata))))
    ## state.names <- rownames(z$trans)
    z
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

## Collapse across states
## issue with variance calculations...
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
    Pdisc <- sum4(object$Pdisc,index)
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
    z <- object
    z$P <- P
    z$L <- L
    z$Pu <- Pu
    z$Lu <- Lu
    z$Pdisc <- Pdisc
    z$trans <- trans
    ## costs??
    z$call <- call
    z
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
predict.stratifiedModel <- function (object, type, newdata, ...) {
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

vsplinefun <- function(x,y,...) {
    stopifnot(is.matrix(y))
    splinefuns <- lapply(1:ncol(y), function(j) stats::splinefun(x,y[,j],...))
    function(x) sapply(splinefuns, function(obj) obj(x))
}
makeSplineFun <- function(object,tm,newdata=data.frame(one=1),tmvar="", min.tm=1e-8, log=FALSE) {
    tm <- pmax(tm,min.tm)
    trans <- if (log) base::log else base::identity
    ## itrans <- if (log) base::exp else base::identity
    newdata <- as.data.frame(newdata)
    times <- data.frame(t=tm); if (tmvar != "") names(times) <- tmvar
    hazard <- lapply(1:nrow(newdata), function(i)
        splinefun(trans(tm),
                  predict(object,type="hazard",newdata=merge(newdata[i,,drop=FALSE], times))))
    gradh <- lapply(1:nrow(newdata), function(i)
        vsplinefun(trans(tm),
                   predict(object,type="gradh",newdata=merge(newdata[i,,drop=FALSE], times))))
    out <- list(
        object=object,
        hazard=function(t) sapply(hazard, function(obj) obj(trans(t))),
        gradh=function(t) t(sapply(gradh, function(obj) obj(trans(t)))))
    class(out) <- "SplineFun"
    out
}
predict.SplineFun <- function(object, newdata, type=c("hazard","gradh"), tmvar="", ...) {
    type <- match.arg(type)
    time <- newdata[[tmvar]]
    stopifnot(all(time[1] == time[-1]))
    if (type=="hazard") object$hazard(time[1]) else object$gradh(time[1])
}
vcov.SplineFun <- function(object) vcov(object$object)
coef.SplineFun <- function(object) coef(object$object)

## Non-parametric baseline: SDE approach due to Ryalen and colleagues
markov_sde <- function(models, trans, newdata, init=NULL, nLebesgue=1e4+1, los=FALSE, nOut=300,
                        weights=1) {
    transfun <- function(tmat) {
        indices <- sort(as.vector(tmat)); indices <- setdiff(indices,NA)
        nStates <- nrow(tmat)
        out <- do.call(rbind,
                lapply(indices, function(i) {
                    index2 <- which(tmat == i)
                    from <- (index2-1) %% nStates +1
                    to <- (index2-1) %/% nStates + 1
                    data.frame(from=from,to=to)
                }))
        matrix(as.integer(as.matrix(out)),nrow(out))-1L
    }
    ## TODO check parameters
    nStates <- nrow(trans)
    nTrans <- sum(!is.na(trans))
    if (is.null(init)) {
        init <- c(1,rep(0,nStates-1))
    }
    stopifnot(length(init)==nStates)
    ## init <- rep(init,nrow(newdata))
    n <- sum(sapply(models,attr,"orig.max.clust"))
    cumHazList <- lapply(models, function(object)
        -t(log(predict(object, newdata=newdata, se=FALSE)$S0)))
    timesList <- lapply(models, function(model) model$cum[,1])
    eventTimes <- times <- sort(unique(unlist(timesList)))
    hazList <- lapply(1:nrow(newdata), function(i) {
        hazMatrix <- matrix(0,nrow=nTrans,ncol=length(times))
        for (j in 1:nTrans) {
            hazMatrix[j,match(timesList[[j]],times)] <- diff(c(0,cumHazList[[j]][,i]))
        }
        hazMatrix
    })
    hazMatrix <- do.call(rbind,hazList)
    if (length(weights)==1 && weights != 0)
        weights <- rep(weights,nrow(newdata))/(weights*nrow(newdata))
    vcov <- matrix(1,1,1)
    if (!los) {
        out <- .Call("plugin_P_by", n, nrow(newdata), hazMatrix, init, transfun(trans), weights,
                     vcov,
                     PACKAGE="rstpm2")
        out$P <- out$X
        out$P.se <- sqrt(out$variance)
    } else {
        out <- .Call("plugin_P_L_by",
                     n, nrow(newdata), hazMatrix, init, transfun(trans), times, weights, nOut,
                     vcov, nLebesgue,
                     PACKAGE="rstpm2")
        PIndex <- 1:(nrow(out$X)/2)
        tr <- function(x) array(as.vector(x), dim=c(nStates,nrow(newdata),ncol(x)))
        out$P <- tr(out$X[PIndex,])
        out$L <- tr(out$X[-PIndex,])
        out$P.se <- sqrt(tr(out$variance[PIndex,]))
        out$L.se <- sqrt(tr(out$variance[-PIndex,]))
    }
    out$n <- n
    out$times <- if (los) out$time else times
    out$newdata <- newdata
    out$trans <- trans
    out$los <- los
    out$init <- init
    out$weights <- weights
    class(out) <- "markov_sde"
    if (!all(weights==0)) {
        stand <- out # copy! This may be a bad idea...
        stand$newdata <- newdata
        if (!los) {
            stand$P <- stand$Y
            stand$P.se <- sqrt(stand$varY)
        } else {
            PIndex <- 1:(nrow(out$Y)/2)
            stand$P <- stand$Y[PIndex,]
            stand$L <- stand$Y[-PIndex,]
            stand$P.se <- sqrt(stand$varY[PIndex,])
            stand$L.se <- sqrt(stand$varY[-PIndex,])
        }
        ## tidy up
        stand$X <- stand$Y
        stand$variance <- stand$varY
        out$X <- out$variance <- out$Y <- out$varY <- stand$Y <- stand$varY <- NULL
        stand$newdata <- stand$newdata[1,,drop=FALSE]
        class(stand) <- "markov_sde"
        out$stand <- stand
    }
    out
}
standardise.markov_sde <- function(x, ...) {
    x$stand
}

plot.markov_sde <- function(x, y, stacked=TRUE, which=c("P","L"), index=NULL,
                            xlab="Time", ylab=NULL, col=2:6, border=col,
                            ggplot2=FALSE, lattice=FALSE, alpha=0.2,
                            strata=NULL,
                            ...) {
    stopifnot(inherits(x,"markov_sde"))
    which <- match.arg(which)
    if (!missing(y)) warning("y argument is ignored")
    ## ylab defaults
    if (is.null(ylab))
        ylab <- if(which=='P') "Probability" else "Length of stay"
    if (ggplot2)
        ggplot.markov_msm(x, which=which, stacked=stacked, xlab=xlab, ylab=ylab,
                          alpha=alpha, ...)
    else if (lattice)
        xyplot.markov_msm(x, which=which, stacked=stacked, xlab=xlab, ylab=ylab,
                          col=col, border=border, strata=strata, ...)
    else {
        ## Is the next statement correct??
        if (is.null(index) && nrow(x$newdata)>1) {
            warning("More than one set of covariates; defaults to weighted estimator")
            x <- x$stand # Warning: replacement
            index <- 1
        }
        if (is.null(index)) index <- 1
        df <- merge(x$newdata[index,,drop=FALSE], as.data.frame(x))
        states <- unique(df$state)
        if (stacked) {
            out <- graphics::plot(range(x$times, na.rm=TRUE),0:1, type="n", xlab=xlab, ylab=ylab, ...)
            lower <- 0
            for (i in length(states):1) { # put the last state at the bottom
                df2 <- df[df$state==states[i],]
                if (length(lower)==1) lower <- rep(0,nrow(df2))
                upper <- lower+df2[[which]]
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

as.data.frame.markov_sde <- function(x, row.names=NULL, optional=NULL, ci=TRUE,
                                      P.conf.type="logit", L.conf.type="log",
                                      P.range=c(0,1), L.range=c(0,Inf),
                                      ...) {
    if (any(x$weights<0)) {
        P.conf.type <- L.conf.type <- "plain"
        P.range <- L.range <- c(-Inf,Inf)
    }
    .id. <- 1:nrow(x$newdata)
    nStates <- nrow(x$trans)
    state.names <- rownames(x$trans)
    stateNames <- if (!is.null(rownames(x$trans))) rownames(x$trans) else 1:nrow(x$trans)
    out <- expand.grid(state=stateNames, .id.=.id., time=x$times)
    out <- cbind(x$newdata[out$.id.,],out)
    names(out)[1:ncol(x$newdata)] <- colnames(x$newdata)
    out$P <- as.vector(x$P)
    out$P.se <- as.vector(x$P.se)
    if (ci) {
        tmp <- surv.confint(out$P,out$P.se, conf.type=P.conf.type, min.value=P.range[1], max.value=P.range[2])
        out$P.lower <- tmp$lower
        out$P.upper <- tmp$upper
    }
    if (x$los) {
        out$L <- as.vector(x$L)
        out$L.se <- as.vector(x$L.se)
        if (ci) {
            tmp <- surv.confint(out$L,out$L.se, conf.type=L.conf.type, min.value=L.range[1], max.value=L.range[2])
            out$L.lower <- tmp$lower
            out$L.upper <- tmp$upper
        }
    }
    out <- out[order(out$.id.,out$state,out$time),]
    out$.id. <- NULL
    if(!is.null(row.names)) rownames(out) <- row.names
    out
}

## Parametric bootstrap
pboot = function(object, m, parallel=FALSE, mc.cores=parallel::detectCores(), ...) {
    ## do this once
    newx = object$x
    coefs <- lapply(object$x, function(xi) coef(xi)) # note this trick!
    vcovs <- lapply(object$x, function(xi) vcov(xi)) # note this trick!
    ## vcov = lapply(object$x, vcov) # fails:(
    ## do the following m times
    parallel::mclapply(1:m , function(i) {
        ## update the coefs
        for (j in 1:length(newx))
            coef(newx[[j]]) = mvtnorm::rmvnorm(1,coefs[[j]],vcovs[[j]])
        ## update the object
        update(object,x=newx)
    }, mc.cores=mc.cores)
}
