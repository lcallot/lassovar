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
	
	dd <- dim(as.matrix(fc.data))

	if(1%in%dd){
		fc.fit	<-cbind(1,matrix(fc.data,nrow=1))%*%coef(object)
		#names(fc.fit)<-object$var.names
	}	
	if(!(1%in%dd)){
		fc.fit	<-cbind(1,fc.data)%*%coef(object)
	}
		colnames(fc.fit)<-object$var.names
	return(fc.fit)
}
