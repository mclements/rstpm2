setClass("tvcCoxph", contains="mle2")
## setClass("tvcCoxphList", contains="list")

cox.tvc <- function(obj,var=NULL,method="logt") {
    stopifnot(attr(obj$y,"type") == "right")
    if (is.null(var))
        return(as(lapply(attr(obj$terms,"term.labels"),
                      function(name) cox.tvc(obj, var=name,method=method)),
                  "tvcCoxList"))
    method <- match.arg(method)
    y <- as.matrix(obj$y)
    time <- y[,1]
    status <- y[,2]
    X <- as.matrix(model.matrix(obj))
    index <- order(time,-status)
    X <- X[index, , drop = FALSE]
    time <- time[index]
    status <- status[index]
    k <- match(var,attr(obj$terms,"term.labels"))
    beta <- c(coef(obj),0)
    names(beta) <- c(names(coef(obj)),sprintf("%s:log(t)",var))
    minuslogl <- function(beta) -.Call("test_cox_tvc2",
                                   list(time=time,event=status,X=X,beta=beta,k=k-1),package="rstpm2")
    gr <- function(beta) -.Call("test_cox_tvc2_grad",
                                list(time=time,event=status,X=X,beta=beta,k=k-1),package="rstpm2")
    parnames(minuslogl) <- parnames(gr) <- names(beta)
    fit <- mle2(start=beta,
                 minuslogl = minuslogl,
                 gr = gr,
                 method="BFGS")
    fit@data <- list(object=obj,k=k,var=var)
    ## return the mle2 object
    as(fit,"tvcCoxph")
}

setMethod("show", "tvcCoxph", function(object) print(summary(object)))

setMethod("plot", signature(x="tvcCoxph", y="missing"),
          function(x,y,
                   add=FALSE,rug=!add,
                   type="l",lty=c(1,2,2),col=1,ylab="Effect",log="y",...) {
              obj <- x@data$object
              y <- as.matrix(obj$y)
              time <- y[,1]
              status <- y[,2]
              timek <- seq(min(time[status==1]), max(time[status==1]), length=301)
              Xk <- cbind(1,log(timek))
              k <- x@data$k
              index <- c(k,length(coef(x)))
              betak <- coef(x)[index]
              fitted <- as.vector(Xk %*% betak)
              gd <- t(Xk)
              se.fit <- sqrt(colSums(gd * (vcov(x)[index,index] %*% gd)))
              if (add) {
                  matlines(timek,exp(fitted+cbind(0,-1.96*se.fit,1.96*se.fit)),lty=lty,col=col,...)
              } else {
                  matplot(timek,exp(fitted+cbind(0,-1.96*se.fit,1.96*se.fit)),type=type,lty=lty,col=col,ylab=ylab,log=log,...)
              }
              abline(h=exp(coef(obj)[k]),lty=3)
              if (rug)
                  rug(time[status==1])
          })
