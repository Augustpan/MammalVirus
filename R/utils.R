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