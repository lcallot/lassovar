summary.lassovar<-function(object,...)
{
	if(!(class(object)=='lassovar'))stop('Object is not of the lassovar class')

	argList	<- list(...)
	ic	<-argList$ic
	short	<-argList$short

	if(is.null(ic))ic<-'BIC'
	if(is.null(short))short<-FALSE
	
	ic.mod	<-object[[ic]]
	
	
	cat('Call:\n',sep='')
	print(object$call)
	cat('\n')

	cat('Model estimated equation by equation','\n',sep='')
	cat('Selection criterion: ',ic,'\n',sep='')				
	cat('Estimator: ',object$estimator,'\n',sep='')
	if(!is.null(object$ada.w)){cat('Adaptive weights Estimator: ',object$ada.w$ada,'\n',sep='')}


	T	<-nrow(object$y)
	N	<-ncol(object$y)
	cat('Dimensions: T = ',T,'  N = ',N,'\n',sep='')
	nbr.nz	<-sum(colSums(ic.mod$coefficients[-1,]!=0))
	if(object$estimator!='ols')
		{
		cat('\nTotal number of variables selected:',nbr.nz,' (',100*nbr.nz/length(ic.mod$coef[-1,]),'% of candidates)\n',sep='')
		eq.sum	<-cbind(ic.mod$lambda,colSums(ic.mod$coefficients[-1,]!=0),ic.mod$rss/T)

		rownames(eq.sum)<-object$var.names
		colnames(eq.sum)<-c('Lambda','non-zero','resid var')
		if(!short){print('Model summary statistics:');print(eq.sum);}
		if(short){cat('Summary statistics','\n','Average for all equations:','\n',sep='');print(colMeans(eq.sum));}
		}

	invisible(eq.sum)
}
