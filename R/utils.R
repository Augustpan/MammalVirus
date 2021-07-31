#' AICc from gam or lm object
#'
#' Helper function compute AICc from a gam or lm object
#' 
#' @param object A model returned by a gam or lm fit
#' 
#' @return an AICc value
#' 
#' @examples
#' # some simulated data with one cov
#' library(mgcv)
#' dat <- gamSim(6,n=20,scale=.5)[,1:2]
#' m <- gam(y~s(x0), data=dat)
#' AICc(m)
AICc <- function(object){
    k <- attributes(logLik(object))$df
    aic <- stats::AIC(object)
    n <- nrow(object$model)
    if(class(object)[1]=="marssMLE") n <- object$samp.size
    return(aic+(2*k^2+2*k)/(n-k-1))
}

####################################
####################################
#### dev_expl()

#' @title Calculate the percent deviance explained by a GLM or GAM
#' @description This function calculates the percent deviance explained by a model, such as a generalised linear model (GLM) or generalised additive model (GAM).
#' @param model A model object.
#' @details Percent deviance explained is given by \code{(model$null.deviance - model$deviance)/model$null.deviance) * 100}.
#' @return The function returns the percent deviance explained as a number.
#' @examples
#' n <- 100
#' x <- stats::runif(n, 0, 10)
#' y <- stats::rpois(n, x*2)
#' m1 <- stats::glm(y ~ x)
#' dev_expl(m1)
#' m2 <- mgcv::gam(y ~ s(x))
#' dev_expl(m2)
#'
#' @author Edward Lavender
#' @export
#'

dev_expl <- function(model){
    if(is.null(model$null.deviance) | is.null(model$deviance)){
        stop("model$null.deviance or model$deviance is NULL.")
    }
    return((model$null.deviance - model$deviance)/model$null.deviance * 100)
}