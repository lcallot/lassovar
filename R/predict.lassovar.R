#' predict.lassovar
#'
#' @param object A lassovar object.
#' @param fc.data A data frame with the new data.
#' @param ... unused
#'
#' @method predict lassovar 
#' @export
predict.lassovar <-
function(object,fc.data,...)
{
	fc.fit	<-cbind(1,fc.data)%*%coef(object)
	names(fc.fit)<-object$var.names
	
	return(fc.fit)
}
