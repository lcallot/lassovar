#' Estimates a Vector Autoregressive model by mean of Lasso or adaptive Lasso
#'
#' @param dat a data.frame containing the series.
#' @param lags the number of consecutive lags in the model. 
#'
#' @param ic Optional, the information criterion to use for selecting the penalty parameter. Default to BIC.
#' @param ada Optional, initial estimator for the adaptive Lasso. Can be chosen between ols, lasso, ridge. 
#' @param rest Optional, prespecified restrictions on the design matrix. Can take values in 'none', 'ar' for autoregressive models, and 'covrest' for covariance-style restrictions. 
#' @param exo Optional, a data frame of same length as dat containing exogenous variables.   
#' @param mc Optional, default to FALSE. Should the equation-by-equation estimation be parallelized using multicore?
#' @param ncores Optional, the number of cores to use in case of parallelization.
#' @param dfmax Optional, the maximum number of variables in the model excluding the intercept. An option of glmnet, it exits the algorithm when the penalty is small enough that more than dmax variables are included in the model. Incresease the speed tremendously for large VARs.
#'
#'
#' @return A lassovar object.
#' 
#' @examples
#' \dontrun{
#' dat <- data.frame(matrix(rnorm(1000,ncol=10)))
#' lv.mod <- lassovar(dat)
#' }
#'
#' @export
lassovar<-function(dat,lags=1,...)
{
	argList	<- list(...)
	
	if(!(is.double(lags)&lags>0))stop('something wrong with the lag argument, an integer greater than zero would be appreciated')

	# Checking if dat is a data frame or coercing to it.
	if(!is.data.frame(dat)) dat<-as.data.frame(dat)
	#Checking if exo exists and if so coercing to dataframe.
	if(!is.null(argList$exo))exo<-as.data.frame(argList$exo) else exo <-NULL	

	if(!is.null(argList$post))post.ols<-argList$post else post.ols <-FALSE
	
	
	if(!is.null(argList$ic))ic<-argList$ic			else ic	 <- 'BIC'
	# multicore ?
	if(!is.null(argList$mc))mc<-argList$mc			else mc	<-FALSE
	if(!is.null(argList$ncores))ncores<-argList$ncores	else ncores	<-NULL
	# restrictions ?
	if(!is.null(argList$rest))rest<-argList$rest	else rest	<-'none'
	# ??	
	if(!is.null(argList$dfmax))dfmax<-argList$dfmax	else dfmax	<-ncol(dat)+1

	y.var	<-.mkvar(dat,lags=lags,horizon=1,exo=exo)	
	
	if(!is.null(argList$adaptive)){
		ada<-argList$adaptive
		cat('initial estimator for the adapive lasso: ',ada,'\n',sep='')
		if(!(ada%in%c('ols','group','lasso','ridge'))){stop('unsupported first step procedure')}
		else{
			for(a in ada)
				{
			if(a=='lasso')	ada.w<-.ada.las.weights(y.var$y,y.var$x,a,ic=ic,mc=mc,ncores=ncores,rest=rest,dfmax=dfmax)
			if(a=='ols')	ada.w<-.ada.ols.weights(y.var$y,y.var$x,a,rest=rest,mc=mc,ncores=ncores)
			if(a=='ridge')	ada.w<-.ada.ridge.weights(y.var$y,y.var$x,a,rest=rest,mc=mc,ncores=ncores,dfmax=dfmax)
			if(a=='group')	ada.w<-.ada.grp.weights(y.var$y,y.var$x,a,rest=rest)
				}
		}
	}
	else ada.w<-NULL
	
	#cat('Estimating the Final Lasso: ','\n',sep='')
	las.mod<-.lassovar.eq(y.var$y,y.var$x,ada.w,ic=ic,mc=mc,ncores=ncores,rest=rest,dfmax=dfmax)


	if(post.ols){	las.mod$post <- .post.ols(y.var$y,y.var$x,sel.pars=las.mod$coefficients!=0,mc=mc,ncores=ncores)}		

	las.mod$nbreq	<-ncol(dat)
	las.mod$call	<-match.call()
	return(las.mod)
}









