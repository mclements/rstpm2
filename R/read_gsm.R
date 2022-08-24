#' Utility to find where a formula call matches a name.
#' This may be cleaner than using grep on strings:)
#' @param name quoted name to match
#' @param x right-hand-side of a formula call
#' @return index of the matching positions
grep_call = function(name,x) {
    local_function = function(x)
        if(length(x)==1) x==name else any(sapply(x, local_function))
    which(sapply(x, local_function))
}

#' Extract design information from an stpm2/gsm object and newdata
#' for use in C++
#' @param object stpm2/gsm object
#' @param newdata list or data-frame used for evaluation
#' @param inflate double value to inflate minimum and maximum times for root finding
#' @return list that can be read by `gsm ssim::read_gsm(SEX args)` in C++
#' @rdname gsm_design
#' @importFrom stats predict
#' @export
gsm_design = function(object, newdata, inflate=100) {
    stopifnot(inherits(object, "stpm2"),
              is.list(newdata),
              is.numeric(inflate),
              length(inflate) == 1)
    ## Assumed pattern:
    ## (ns|nsx)(log(timeVar),knots,Boundary.knots,centre=FALSE,derivs=(c(2,2)|c(2,1)))
    terms = attr(object@model.frame, "terms")
    index = grep_call(object@timeVar, attr(terms, "variables"))
    if(length(index)==0) stop("No timeVar in the formula")
    if(length(index)>1)
        stop("Current implementation does not allow for more than one timeVar expression")
    mycall = attr(terms, "predvars")[[index]]
    df = length(c(mycall$knots, mycall$Boundary.knots)) - 1
    stopifnot(mycall[[1]] == quote(nsx) || mycall[[1]] == quote(ns),
              length(mycall[[2]])>1,
              mycall[[2]][[1]] == quote(log), # assumes log
              mycall[[2]][[2]] == object@timeVar, # what about a scalar product or divisor?
              is.null(mycall$deriv) || (mycall$derivs[1] == 2 && mycall$derivs[2] %in% 1:2),
              mycall$centre == FALSE) # doesn't allow for centering
    cure = !is.null(mycall$derives) && all(mycall$derivs == c(2,1))
    time = object@args$time
    newdata[[object@timeVar]] = mean(time) # NB: time not used
    Xp = predict(object, newdata=newdata, type="lpmatrix")
    index2 = 1:(ncol(Xp)-df)
    etap = drop(Xp[, index2, drop=FALSE] %*% coef(object)[index2])
    q_const = attr(nsx(log(mean(time)), knots=mycall$knots, Boundary.knots=mycall$Boundary.knots,
                       intercept=mycall$intercept),
                   "q.const")
    list(type="gsm",
         link_name=object@args$link,
         call = mycall,
         tmin = min(time), # not currently used?
         tmax = max(time),
         etap=etap, 
         knots=mycall$knots,
         Boundary_knots=mycall$Boundary.knots,
         intercept=as.integer(mycall$intercept),
         gamma=coef(object)[-index2],
         q_const = q_const,
         cure = as.integer(cure),
         log_time=TRUE,
         inflate=as.double(inflate))
}
