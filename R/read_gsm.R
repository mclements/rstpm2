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
    ## Assumed patterns:
    ## timeEffect := (ns|nsx)(log(timeVar),knots,Boundary.knots,centre=FALSE,derivs=(c(2,2)|c(2,1)))
    ## effect := timeEffect | otherEffect:timeEffect | timeEffect:otherEffect
    terms = attr(object@model.frame, "terms")
    factors = attr(terms, "factors")[-1,,drop=FALSE]
    variables = attr(terms, "variables")[-(1:2)]
    predvars = attr(terms, "predvars")[-(1:2)]
    indices = grep_call(object@timeVar, variables)
    if(length(indices)==0) stop("No timeVar in the formula -- unexpected error")
    index_time_variables = grep_call(object@timeVar, variables) # time variables
    index_time_effects = grep(object@timeVar,colnames(factors)) # components in the rhs with time variables
    ## We need to know how wide is each term
    nms = names(coef(object))
    term.labels = attr(terms, "term.labels")
    coef_index <-
        sapply(strsplit(nms, ":"), function(c) {
            if (length(c)>2)
                stop("current implementation only allows for main effects and two-way interactions")
            pmatchp = function(x, table)
                !is.na(pmatch(x, table))
            if(length(c)==1) {
                for (i in seq_along(term.labels)) {
                    t = term.labels[i]
                    if (pmatchp(t,c))
                        return(i)
                }
                if (c == "(Intercept)")
                    return(0)
                return(-1)
            }
            ## => length(c) == 2
            term.labels.split = strsplit(term.labels, ":")
            for (i in seq_along(term.labels)) {
                t = term.labels.split[[i]]
                if (all(pmatchp(t,c)))
                    return(i)
            }
            return(-1)
        })
    parse_ns = function(mycall,x,index_time_effect) {
        df = length(c(mycall$knots, mycall$Boundary.knots)) - 1
        stopifnot(mycall[[1]] == quote(nsx) || mycall[[1]] == quote(ns),
                  length(mycall[[2]])>1,
                  mycall[[2]][[1]] == quote(log), # assumes log
                  mycall[[2]][[2]] == object@timeVar, # what about a scalar product or divisor?
                  is.null(mycall$deriv) || (mycall$derivs[1] == 2 && mycall$derivs[2] %in% 1:2),
                  mycall$centre == FALSE) # doesn't allow for centering
        cure = !is.null(mycall$derives) && all(mycall$derivs == c(2,1))
        time = object@args$time
        q_const = attr(nsx(log(mean(time)), knots=mycall$knots,
                           Boundary.knots=mycall$Boundary.knots,
                           intercept=mycall$intercept),
                       "q.const")
        list(call = mycall,
             knots=mycall$knots,
             Boundary_knots=mycall$Boundary.knots,
             intercept=as.integer(mycall$intercept),
             gamma=coef(object)[which(coef_index %in% index_time_effect)],
             q_const = q_const,
             cure = as.integer(cure),
             x=x)
    }
    time = object@args$time
    newdata[[object@timeVar]] = mean(time) # NB: time not used
    Xp = predict(object, newdata=newdata, type="lpmatrix")
    index2 = which(!(coef_index %in% index_time_effects))
    etap = drop(Xp[, index2, drop=FALSE] %*% coef(object)[index2])
    list(type="gsm",
         link_name=object@args$link,
         tmin = min(time), # not currently used?
         tmax = max(time),
         inflate=as.double(inflate),
         etap=etap,
         coefp = coef(object)[index2], # for debugging
         log_time=TRUE,
         terms =
             lapply(index_time_effects,
                    function(i) {
                        j = which(factors[,i] != 0)
                        if (length(j)==1)
                            return(parse_ns(predvars[[j]], rep(1, nrow(newdata)), i))
                        else {
                            if(length(j)>3)
                                stop("Current implementation only allows for two-way interaction terms")
                            if (j[1] %in% index_time_variables)
                                return(parse_ns(predvars[[j[1]]], eval(predvars[[j[2]]], newdata), i))
                            else return(parse_ns(predvars[[j[2]]], eval(predvars[[j[1]]], newdata), i))
                        }
                    })
         )
}
