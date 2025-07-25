.onLoad <- function(libname, pkgname) {
    ## options(warn = 3, error = utils::recover)
}

.onUnload <- function (libpath) {
  library.dynam.unload("rstpm2", libpath)
}
