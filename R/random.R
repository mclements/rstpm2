## Drop-in replacement functions (not exported)
rgompertz = function(n, shape, rate) {
    u <- 1-runif(n)
    i <- shape >= 0 | u>=exp(rate/shape)
    y <- rep(Inf, n)
    if (any(i))
        y[i] <- log(1 - shape*log(u[i])/rate) / shape
    y
}
rllogis = function(n, shape, scale) {
    u <- runif(n)
    scale*(u/(1-u))^(1/shape)
}
