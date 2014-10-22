tvc.coxph <- function(obj,var,method="logt") {
    stopifnot(attr(obj$y,"type") == "right")
    stopifnot(method == "logt")
    y <- as.matrix(obj$y)
    time <- y[,1]
    status <- y[,2]
    X <- as.matrix(model.matrix(obj))
    index <- order(time,-status)
    X <- X[index, , drop = FALSE]
    time <- time[index]
    status <- status[index]
    k <- match(attr(obj$terms,"term.labels"),var)
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
                 method="BFGS", hessian=TRUE)
    fit@data <- list(object=obj)
    ## plot the result
    betak <- beta*0
    index2 <- c(k,length(betak))
    betak[index2] <- coef(fit)[index2]
    Xk <- cbind(X,log(time))
    Xk[,-index2] <- 0
    Xk[,k] <- 1
    fitted <- as.vector(Xk %*% betak)
    gd <- t(Xk)
    se.fit <- sqrt(colSums(gd * (vcov(fit) %*% gd)))
    matplot(time,exp(fitted+cbind(0,-1.96*se.fit,1.96*se.fit)),type="l",lty=c(1,2,2),col=1,ylab="Effect",log="y")
    abline(h=exp(coef(obj)[k]),lty=3)
    fit
}
