#' The residuals method of class lassovar.
#' 
#' @aliases residuals resid
#' 
#' @param object A lassovar object.
#' @param ... use short=TRUE for a shorter printout.
#' @method residuals lassovar 
#' @export
residuals.lassovar <- function(object,...){
	if(!(class(object)=='lassovar'))stop('Object is not of the lassovar class')

	res <- object$y-cbind(1,object$x)%*%object$coefficients

	return(res)
}
