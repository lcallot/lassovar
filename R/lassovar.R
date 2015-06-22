#' Fits a Vector Autoregressive model by l1 penalization.
#' 
#' @description Fits a Vector Autoregressive model by mean of Lasso or adaptive Lasso. The parameters are estimated using \pkg{glmnet}, the penalty parameter is selected using information criteria. The VAR is constructed from the \code{data.frame} provided, exogenous variables and various deterministic specifications are possible in options. 
#'
#' @param dat a data.frame containing the series.
#' @param exo Optional, a data frame of same length as dat containing exogenous variables. The exogenous variables are not lagged.   
#' @param lags the number of consecutive lags in the model. 
#'
#' @param ic Optional, the information criterion to use for selecting the penalty parameter. Default to BIC.
#' @param adaptive Optional, initial estimator for the adaptive Lasso. Can be chosen between: none (default, plain Lasso), ols, lasso, ridge. 
#' @param mc Optional, default to FALSE. Should the equation-by-equation estimation be parallelized using mclapply from the parallel package?
#' @param ncores Optional, the number of cores to use in case of parallelization.
#' @param dfmax Optional, the maximum number of variables in the model excluding the intercept. An option of glmnet, it exits the algorithm when the penalty is small enough that more than dmax variables are included in the model. Incresease the speed tremendously for large VARs.
#' @param post Optional, Should a post Lasso OLS be estimated, default FALSE.
#' @param horizon Estimate a h-step ahead VAR, useful for direct forecasting. Default = 1.
#' @param trend Should a linear trend be included in the model. 
#'
#' @return 
#' A list with S3 class \pkg{lassovar}.
#' \describe{
#' 	\item{call}{The call. }
#' 	\item{var.names}{ A vector with the names of the endogenous variables. }
#' 	\item{ada.w}{ A list containing the initial estimator in case of an adaptive Lasso, NULL otherwise.  }
#' 	\item{x}{ A data.frame containing the right hand side variables. }
#' 	\item{y}{ A data.frame containing the left hand side variables. }
#' 	\item{coefficients}{  The parameters of the VAR. }
#' 	\item{RSS}{ The residual sum of squares.  }
#' 	\item{Lambda}{ A vector with the penalty parameter selected for each equation. }
#' 	\item{spectes}{ A matrix containing the results of specification tests (TO BE COMPLETED & ADD SPECTESTS TO SUMMARY). }
#' 	\item{estimator}{ A string containing the name of the estimator (lasso or adaptive lasso). }  
#' 	\item{ic}{ A string with the name of the information criterion used to select the penlaty parameters. } 
#' 	\item{nbreq}{ The number of equations. }
#' }
#' 
#' @examples
#' \dontrun{
#' dat <- data.frame(matrix(rnorm(100),ncol=5))
#' lv.mod <- lassovar(dat,lags=1)
#' lv.mod.ada <- lassovar(dat,lags=2,adaptive ='ols')
#' }
#'
#' @export
lassovar<-function(dat,exo=NULL,lags=1,ic=c('BIC','AIC'),adaptive=c('none','ols','lasso','group','ridge'),post=FALSE,mc=FALSE,ncores=NULL,dfmax=NULL,horizon=1,trend=FALSE)
{
	
	# matching the multiple choice arguments. 
	ic<-match.arg(ic)
	adaptive <- match.arg(adaptive)
	
	if(!((lags==as.integer(lags))&lags>0))stop('Invalid \'lags\' argument. Positive integer required.')

	# Checking if dat is a data frame or coercing to it.
	if(!is.data.frame(dat)) dat<-as.data.frame(dat)
	
	#Checking if exo exists and if so coercing to dataframe.
	if(!is.null(exo))exo<-as.data.frame(exo) else exo <-NULL	
	
	#if(!is.null(post))post<-post else post <-FALSE
	#if(!is.null(mc))mc<-mc else mc	<-FALSE
	
	# multicore ?
	if(!is.null(ncores))ncores<-ncores	else ncores	<-NULL
	
	# maxdegrees of freedom ?	
	if(!is.null(dfmax))dfmax<-as.integer(dfmax)	else dfmax	<-ncol(dat)*lags

	y.var	<-.mkvar(dat,lags=lags,horizon=1,exo=exo,trend=trend)	
	
	if(adaptive!='none'){
		cat('initial estimator for the adapive lasso: ',adaptive,'\n',sep='')
		for(a in adaptive){
			if(a=='lasso')	ada.w<-.ada.las.weights(y.var$y,y.var$x,a,ic=ic,mc=mc,ncores=ncores,dfmax=dfmax,trend=trend)
			if(a=='ols')	ada.w<-.ada.ols.weights(y.var$y,y.var$x,a,mc=mc,ncores=ncores)
			if(a=='ridge')	ada.w<-.ada.ridge.weights(y.var$y,y.var$x,a,mc=mc,ncores=ncores,dfmax=dfmax,trend=trend)
			if(a=='group')	ada.w<-.ada.grp.weights(y.var$y,y.var$x,a,trend)
		}
	}
	else ada.w<-NULL
	
	#cat('Estimating the Final Lasso: ','\n',sep='')
	las.mod<-.lassovar.eq(y.var$y,y.var$x,ada.w,ic=ic,mc=mc,ncores=ncores,dfmax=dfmax,trend=trend)


	if(post){	las.mod$post <- .post.ols(y.var$y,y.var$x,sel.pars=las.mod$coefficients!=0,mc=mc,ncores=ncores)}		

	las.mod$trend 	<-trend
	las.mod$ic  	<-ic
	las.mod$nbreq	<-ncol(dat)
	las.mod$call	<-match.call()
	return(las.mod)
}









