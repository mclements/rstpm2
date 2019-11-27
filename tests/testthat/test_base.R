library(rstpm2)

## for coping with weird test behaviour from CRAN and R-devel
.CRAN <- TRUE
## pstpm2+frailty models are slow
slow <- FALSE

expect_eps <- function(expr, value, eps=1e-7)
    expect_lt(max(abs(expr-value)),eps)

context("stpm2")
##
test_that("base", {
    fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer)
    expect_eps(coef(fit)[2], -0.361403, 1e-5)
})
test_that("NelderMead", {
    fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,
                 control=list(optimiser="NelderMead"))
    expect_eps(coef(fit)[2], -0.3608321, 1e-5)
})
test_that("tvc", {
    ## main effect IS needed for the parametric case (else constrained to be 0 at first event time)
    fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,tvc=list(hormon=2))
    expect_eps(coef(fit)[2], -0.9930585, 1e-5)
})
test_that("Comparison with Stata", {
    fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,
                 smooth.formula=~nsx(log(rectime),df=3,stata=TRUE))
    expect_eps(coef(fit)[2], -0.3614357, 1e-5)
})
test_that("Cure", {
    fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,
                 smooth.formula=~nsx(log(rectime),df=3,cure=TRUE))
    expect_eps(coef(fit)[2], -0.3564268, 1e-5)
    fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,cure=TRUE)
    expect_eps(coef(fit)[2], -0.3564268, 1e-5)
})
test_that("probit", {
    fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,link="probit")
    expect_eps(coef(fit)[2], -0.282318, 1e-5)
})


test_that("interval", {
    ## ICPHREG example
    read.textConnection <- function(text, ...) {
        conn <-  textConnection(text)
        on.exit(close(conn))
        read.table(conn, ...)
    }
    hiv <- read.textConnection("0 16 0 0 0 1
15 26 0 0 0 1
12 26 0 0 0 1
17 26 0 0 0 1
13 26 0 0 0 1
0 24 0 0 1 0
6 26 0 1 1 0
0 15 0 1 1 0
14 26 0 1 1 0
12 26 0 1 1 0
13 26 0 1 0 1
12 26 0 1 1 0
12 26 0 1 1 0
0 18 0 1 0 1
0 14 0 1 0 1
0 17 0 1 1 0
0 15 0 1 1 0
3 26 1 0 0 1
4 26 1 0 0 1
1 11 1 0 0 1
13 19 1 0 0 1
0 6 1 0 0 1
0 11 1 1 0 0
6 26 1 1 0 0
0 6 1 1 0 0
2 12 1 1 0 0
1 17 1 1 1 0
0 14 1 1 0 0
0 25 1 1 0 1
2 11 1 1 0 0
0 14 1 1 0 0")
    names(hiv) <- c("Left","Right","Stage","Dose","CdLow","CdHigh")
    hiv <- transform(hiv,Event = ifelse(Left==0,2,
                                 ifelse(Right>=26,0,
                                        3)))
    hiv2 <- transform(hiv,
                      Left = ifelse(Event==2,Right,
                             ifelse(Event==0,Left,
                                    Left)),
                      Right = ifelse(Event==2,Inf,
                              ifelse(Event==0,Inf,
                                     Right)))
    fit <- stpm2(Surv(Left,Right,Event,type="interval")~Stage, data=hiv2, df=2)
    expect_eps(coef(fit)[2], 1.917699, 1e-4)
})


context("stpm2 + frailty")
##
test_that("base", {
    set.seed(12345)
    id <- rep(1:200,each=2)
    x <- rnorm(400)
    mu.cluster <- rnorm(200)
    y <- rexp(400, exp(mu.cluster[id]+x))
    d <- data.frame(y,e=1,x,id)
    fit <- stpm2(Surv(y,e)~x+cluster(id),data=d)
    expect_eps(coef(fit)[2], 0.9256663, 1e-5)
    expect_eps(coef(fit)[6], -0.4698344, 1e-5)
    fit <- stpm2(Surv(y,e)~x+cluster(id),data=d,RandDist="LogN")
    expect_eps(coef(fit)[2], 0.93819784, 1e-5)
    expect_eps(coef(fit)[6], 0.02906888, 1e-5)
})
##
test_that("providing initial values with frailty works", {
  dat <- data.frame(
    y     = c(1, 2, 1, 2, 1, 1), 
    event = c(0, 1, 1, 0, 1, 1),
    grp   = c(1, 1, 2, 2, 3, 3))
  
  fit     <- gsm(Surv(y, event) ~ 1, dat, df = 1, cluster = dat$grp)
  fit_new <- gsm(Surv(y, event) ~ 1, dat, df = 1, cluster = dat$grp, 
                 init = head(coef(fit), -1), logtheta = tail(coef(fit), 1))
  
  expect_equal(coef(fit), coef(fit_new))
})

context("pstpm2")
##
if (!.CRAN)
    test_that("base", {
        fit <- pstpm2(Surv(rectime,censrec==1)~hormon,data=brcancer)
        expect_eps(coef(fit)[2], -0.3650564, 1e-5)
        ## fit <- pstpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,
        ##               control=list(optimiser="NelderMead"))
        ## expect_eps(coef(fit)[2], -0.36581571, 1e-3)
    })

if (slow) {
    test_that("tvc arg", {
        fit <- pstpm2(Surv(rectime,censrec==1)~1,data=brcancer,tvc=list(hormon=-1))
        expect_eps(coef(fit)[1], -1.456845e+00, 1e-5)
        ## main effect is removed by the tvc call
        fit2 <- pstpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,tvc=list(hormon=-1))
        expect_eps(coef(fit), coef(fit2), 1e-10)
        fit3 <- pstpm2(Surv(rectime,censrec==1)~hormon+x1,data=brcancer,tvc=list(hormon=-1))
        expect_eps(coef(fit3)[2], -2.426229e-04, 1e-5)
    })
    test_that("tvc using smooth.formula", {
        ## main effect should be EXCLUDED for the penalised case
        fit <- pstpm2(Surv(rectime,censrec==1)~1,data=brcancer,
                      smooth.formula=~s(log(rectime))+s(log(rectime),by=hormon))
        expect_eps(coef(fit)[1], -1.456845e+00, 1e-5)
    })
    context("pstpm2 + frailty")
    ##
    test_that("base", {
        set.seed(12345)
        id <- rep(1:200,each=2)
        x <- rnorm(400)
        mu.cluster <- rnorm(200)
        y <- rexp(400, exp(mu.cluster[id]+x))
        d <- data.frame(y,e=1,x,id)
        fit <- pstpm2(Surv(y,e)~x+cluster(id),data=d)
        expect_eps(coef(fit)[2], 0.9071705, 1e-5)
        expect_eps(coef(fit)[12], -5.303905e-01, 1e-5)
        fit <- pstpm2(Surv(y,e)~x+cluster(id),data=d,RandDist="LogN")
        expect_eps(coef(fit)[2], 9.334020e-01, 1e-5)
        expect_eps(coef(fit)[6], 1.426793e-02, 1e-5)
    })
}
