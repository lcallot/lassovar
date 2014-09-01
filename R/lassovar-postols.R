
# Post Lasso ols estimator
.post.ols<-
function(y,x,sel.pars,mc,ncores)
{
	post.est<-list('ols'=NULL)

	#ols.b	<-lapply(1:ncol(y),.post.ols.eq,y,x,sel.pars)
	ols.b	<-mclapply(1:ncol(y),.post.ols.eq,y,x,sel.pars,mc.cores=ncores)
	ols.b	<- do.call(cbind,ols.b)

	ols.b <- Matrix(ols.b,sparse=TRUE)
	return(ols.b)
}

# post lasso internal
.post.ols.eq<-	function(i,y,x,sel.pars){
	ols.b 			<- matrix(0,nrow=ncol(x)+1,ncol=1)
	if(sum(sel.pars[-1,i])!=0)ols.b[(1:length(ols.b))*sel.pars[,i]]	<- lm(y[,i]~x[,sel.pars[-1,i]])$coefficient
	else ols.b[1]=mean(y[,i])

	return(ols.b)
}
