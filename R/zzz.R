.onLoad <- function(libname, pkgname) {
    ## trace("formals", tracer = quote({
    ##     cat(">>> formals called in R CMD check with object of class:", class(fun), "\n")
    ##     print(sys.calls())
    ## }), print = FALSE)
    options(warn = 3, error = utils::recover)
}

.onUnload <- function (libpath) {
  library.dynam.unload("rstpm2", libpath)
}
