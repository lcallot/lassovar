#' The summary method of class lassovar.
#' 
#' @param object A lassovar object.
#' @param ... use short=TRUE for a shorter printout.
#' @method summary lassovar 
#' @export
summary.lassovar <- function(object,...){
	if(!(class(object)=='lassovar'))stop('Object is not of the lassovar class')

	argList	<- list(...)
	short	<-argList$short
	if(is.null(short))short<-FALSE
	
	
	cat('Call:\n',sep='')
	print(object$call)
	cat('\n')

	cat('Model estimated equation by equation','\n',sep='')
	cat('Selection criterion: ',object$ic,'\n',sep='')				
	cat('Estimator: ',object$estimator,'\n',sep='')
	if(!is.null(object$ada.w)){cat('Adaptive weights Estimator: ',object$ada.w$ada,'\n',sep='')}
	cat('Deterministics: ',ifelse(object$trend,'interept and trend','intercept'),'.\n',sep='')

	T	<-nrow(object$y)
	N	<-ncol(object$y)
	cat('Dimensions: T = ',T,'  N = ',N,'\n',sep='')
	
	
	if(ncol(object$coefficients)==1) fsums <- sum
	if(ncol(object$coefficients)>1) fsums <- colSums
	
	
	if(!object$trend)nzeq	<-fsums(object$coefficients[-1,]!=0)
	if(object$trend)nzeq	<-fsums(object$coefficients[2:(nrow(object$coefficients)-1),]!=0)
	nbr.nz <- sum(nzeq)
	
	if(object$estimator!='ols')
		{
		cat('\nTotal number of variables selected:',nbr.nz,' (',100*nbr.nz/length(object$coef[-1,]),'% of candidates)\n',sep='')
		eq.sum	<-cbind(matrix(object$lambda,ncol=1),nzeq,matrix(object$RSS/T,ncol=1))

		rownames(eq.sum)<-object$var.names
		colnames(eq.sum)<-c('Lambda','non-zero','resid var')
		if(!short){cat('\n Model summary statistics:\n');print(eq.sum);cat('\nDeterministics not included in the non-zero count.\n')}
		if(short){cat('Summary statistics','\n','Average for all equations:','\n',sep='');print(colMeans(eq.sum));}
		}

	invisible(eq.sum)
}
