#' @importFrom glmnet glmnet
#' @importFrom Matrix Matrix
#' @importFrom parallel mclapply



# The workhorse, fits a lasso possibly adaptive given weights ada.w and selects the best model according both information criteria. 
.lassovar.eq <-
function(y,x,ada.w,degf.type=NULL,ic,mc=FALSE,ncores=1,alpha=1,dfmax)
{
	lasso.eq	<-list('call'=match.call(),'eq'=TRUE,'var.names'=colnames(y),'ada.w'=ada.w,'x'=x,'y'=y,'coefficients'=NULL,'RSS'=NULL,'lambda'=NULL,'spectest'=NULL)	
	all.ic		<-list()

	
	#Estimation with and w/o multicore
	if(!mc){for(i in 1:ncol(y)){ all.ic[[i]]	<-.lv.eq.gn(i,y,x,ada.w,ic=ic,alpha=alpha,dfmax=dfmax)}}
	if(mc){	all.ic<-mclapply(1:ncol(y),.lv.eq.gn,y,x,ada.w,ic=ic,alpha=alpha,dfmax=dfmax,mc.cores=ncores)}


	#Sorting out the IC results
	# This is stoooopid, and costly.
	for(i in 1:ncol(y)){
		lasso.eq$coefficients<-cbind(lasso.eq$coefficients,all.ic[[i]]$coefficients)
		lasso.eq$RSS	<-cbind(lasso.eq$RSS,all.ic[[i]]$rss)
		lasso.eq$lambda	<-cbind(lasso.eq$lambda,all.ic[[i]]$lambda)
		lasso.eq$spectest <-cbind(lasso.eq$spectest,all.ic[[i]]$spectest)
				}
	rm('all.ic')
	gc()


	if(is.null(ada.w))lasso.eq$estimator	<-'Lasso'
	if(!is.null(ada.w)){lasso.eq$estimator	<-'Adaptive Lasso'}

	class(lasso.eq)	<-'lassovar'

return(lasso.eq)
}


# Core function. Estimation of the Lasso and model selection.
.lv.eq.gn<-function(i,y,x,ada.w,ic,alpha,dfmax)
{

	
# Computing exclusion vectors.	
	all.excl<-NULL
	
  # Plain lasso
	if(is.null(ada.w)){
		gn.mod	<-glmnet(x=x,y=y[,i],family='gaussian',exclude=all.excl,alpha=alpha,dfmax=dfmax,standardize=TRUE,type.gaussian='covariance')}
	
	# In case of adaptive lasso
	if(!is.null(ada.w)){
		if(sum(ada.w$w[,i]==Inf)==length(ada.w$w[,i])){gn.mod<-NULL}
		else{
			wzero		<-ada.w$w[,i]
			wzero[which(wzero==Inf)]<-0
			#cat(c(cov.excl,which(ada.w$w[,i]==Inf)),'\n')
			gn.mod	<-glmnet(x=x,y=y[,i],family='gaussian',exclude=c(all.excl,which(ada.w$w[,i]==Inf)),penalty.factor=wzero,alpha=alpha,dfmax=dfmax,standardize=FALSE,type.gaussian='covariance')}
	}
  
	if(!is.null(gn.mod)){lv.ic		<-.ic.modsel(gn.mod,yi=y[,i],x=x,ic=ic,alpha=alpha)}
  
  # Managing the case where the adaptive lasso returned an empty equation
	if(is.null(gn.mod)){
		lv.ic <- list()
		lv.ic$RSS		<- sum((y[,i]-mean(y[,i]))^2)
		lv.ic$coefficients	<- rbind(mean(y[,i]),matrix(0,nrow=ncol(x),ncol=1))
		lv.ic$lambda	<- 0
		lv.ic$spectest 		<- .sptest(matrix(y[,i]-mean(y[,i]),ncol=1),y[,i])
		lv.ic$ic   		<- ic
		lv.ic$ic.all<-NULL
		}

	lv.ic$nbreq	<-ncol(y)

	return(lv.ic)	
}




# After lv is estimated this function computes the Information Criterion and returns the appropriate coefficients, residuals and penalty parameters and post ols (one day).
.ic.modsel<-function(gn.mod,yi,x,ic,alpha)
{	
	fit	<-predict(gn.mod,x)
	res	<-matrix(rep(yi,ncol(fit)),ncol=ncol(fit))-fit

	RSS		<-colSums(res^2)
	if(alpha==1)lv.df	<-gn.mod$df
	if(alpha==0)lv.df	<-.ridge.df(x,gn.mod$lambda)
	
	if(ic=='BIC')ic.all	<- log(RSS/gn.mod$nobs) + lv.df*log(gn.mod$nobs)/gn.mod$nobs
	if(ic=='AIC')ic.all	<- log(RSS/gn.mod$nobs) + lv.df/gn.mod$nobs

	ic.coef	<-rbind(gn.mod$a0,as.matrix(gn.mod$beta))[,which.min(ic.all)]
	ic.rss	<-RSS[which.min(ic.all)]
	ic.lmin	<-gn.mod$lambda[which.min(ic.all)]
	
	# Spectests
	spectest 		<- .sptest(res[,which.min(ic.all)],yi)

	lasso.ic<-list('ic'=ic,'ic.all'=ic.all,'coefficients'=ic.coef,'rss'=ic.rss,'lambda'=ic.lmin,'spectest'=spectest)

return(lasso.ic)
}


.mkvar <-
	function(data.var,lags,horizon,exo=NULL)
	{
		nbrser	<-ncol(data.var)
		
		# The dependent variable
		y	<-data.var
		if(!is.null(colnames(data.var)))varny<-colnames(data.var)
		if(is.null(colnames(data.var)))	varny<-paste('Eq_',1:nbrser,sep='')
		colnames(y)<-varny
		
		# The lags 
		x <- NULL
		varnx <- NULL
		for(l in 1:lags){
			x <- cbind(x,rbind(matrix(NA,ncol=ncol(y),nrow=l+horizon-1),
					  	  as.matrix(head(data.var,-(l+horizon-1)))))
			varnx <- c(varnx,paste(l,'L_',varny,sep=''))
		}
		colnames(x) <- varnx
		
		# Exo variables
		#if(!is.null(exo)){	x<-cbind(x,t(aaply(exo,2,lag,k=1)))}
		if(!is.null(exo)){	x<-cbind(x,as.matrix(exo))}
		
		
		# Trimming  and x
		y <- tail(y,-(lags+horizon-1))
		x <- tail(x,-(lags+horizon-1))
		
		y.var<-list('y'=y,'x'=x,'lags'=lags,'horizon'=horizon,'nbrser'=nbrser)
		return(y.var)
	}



#Takes the matrix of residuals, apply some specification tests (autocorr, heteroskedasticity, normality, ...)
.sptest<-function(mres,yi,lags=min(floor(3*fitdf),nrow(mres)/2),fitdf=1)
	{
	sptest <- list()
	

	BP <- Box.test(mres,lag=lags,type="L",fitdf=floor(fitdf))$p.val
	SW <- shapiro.test(mres)$p.value
	R2 <- 1-((length(mres)*var(mres,na.rm=TRUE))/(var(yi)*length(yi)))

	sptest <- c(BP,SW,R2)

	return(sptest)
	}


