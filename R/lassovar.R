#' Estimates a Vector Autoregressive model by mean of Lasso or adaptive Lasso
#'
#' @param dat a data.frame containing the series.
#' @param exo Optional, a data frame of same length as dat containing exogenous variables.   
#' @param lags the number of consecutive lags in the model. 
#'
#' @param ic Optional, the information criterion to use for selecting the penalty parameter. Default to BIC.
#' @param adaptive Optional, initial estimator for the adaptive Lasso. Can be chosen between: none (default, plain Lasso), ols, lasso, ridge. 
#' @param mc Optional, default to FALSE. Should the equation-by-equation estimation be parallelized using mclapply from the parallel package?
#' @param ncores Optional, the number of cores to use in case of parallelization.
#' @param dfmax Optional, the maximum number of variables in the model excluding the intercept. An option of glmnet, it exits the algorithm when the penalty is small enough that more than dmax variables are included in the model. Incresease the speed tremendously for large VARs.
#' @param post Optional, Should a post Lasso OLS be estimated, default FALSE.
#' @param ... Optional options...
#'
#'
#' @return A lassovar object.
#' 
#' @examples
#' \dontrun{
#' dat <- data.frame(matrix(rnorm(100),ncol=5))
#' lv.mod <- lassovar(dat,lags=1)
#' lv.mod.ada <- lassovar(dat,lags=1,adaptive ='ols')
#' }
#'
#' @export
lassovar<-function(dat,exo=NULL,lags=1,ic=c('BIC','AIC'),adaptive=c('none','ols','lasso','group','ridge'),post=FALSE,mc=FALSE,ncores=NULL,dfmax=NULL,...)
{
	argList	<- list(...)
	
	# matching the multiple choice arguments. 
	ic<-match.arg(ic)
	adaptive <- match.arg(adaptive)
	
	if(!((lags==as.integer(lags))&lags>0))stop('Invalid \'lags\' argument. Positive integer required.')

	# Checking if dat is a data frame or coercing to it.
	if(!is.data.frame(dat)) dat<-as.data.frame(dat)
	
	#Checking if exo exists and if so coercing to dataframe.
	if(!is.null(exo))exo<-as.data.frame(exo) else exo <-NULL	
	
	if(!is.null(post))post<-post else post <-FALSE
	if(!is.null(mc))mc<-mc else mc	<-FALSE
	
	# multicore ?
	if(!is.null(ncores))ncores<-ncores	else ncores	<-NULL
	
	# maxdegrees of freedom ?	
	if(!is.null(dfmax))dfmax<-as.integer(dfmax)	else dfmax	<-ncol(dat)*lags

	y.var	<-.mkvar(dat,lags=lags,horizon=1,exo=exo)	
	
	if(adaptive!='none'){
		cat('initial estimator for the adapive lasso: ',adaptive,'\n',sep='')
		for(a in adaptive){
			if(a=='lasso')	ada.w<-.ada.las.weights(y.var$y,y.var$x,a,ic=ic,mc=mc,ncores=ncores,dfmax=dfmax)
			if(a=='ols')	ada.w<-.ada.ols.weights(y.var$y,y.var$x,a,mc=mc,ncores=ncores)
			if(a=='ridge')	ada.w<-.ada.ridge.weights(y.var$y,y.var$x,a,mc=mc,ncores=ncores,dfmax=dfmax)
			if(a=='group')	ada.w<-.ada.grp.weights(y.var$y,y.var$x,a)
		}
	}
	else ada.w<-NULL
	
	#cat('Estimating the Final Lasso: ','\n',sep='')
	las.mod<-.lassovar.eq(y.var$y,y.var$x,ada.w,ic=ic,mc=mc,ncores=ncores,dfmax=dfmax)


	if(post){	las.mod$post <- .post.ols(y.var$y,y.var$x,sel.pars=las.mod$coefficients!=0,mc=mc,ncores=ncores)}		

	las.mod$ic  	<-ic
	las.mod$nbreq	<-ncol(dat)
	las.mod$call	<-match.call()
	return(las.mod)
}









