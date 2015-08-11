#' Multiple forecasts a VAR using the Lasso or adaptive lasso.
#'
#'
#' @description This function is designed to compute pseudo out-of-sample forecasting experiments. For out-of-sample forecasts use the predict method. The function will produce h-step ahead forecasts either directly, by fitting a h-step ahead VAR, or recursively. The model is trained on fc.train observations, and as many h-step ahead forecasts as possible are computed on the remainder of the sample.   
#' @importFrom parallel mclapply
#'
#' @param dat A data.frame.
#' @param exo Optional, a data frame of same length as dat containing exogenous variables. The exogenous variables are not lagged.   
#' @param fc.train The number of training observations to use.
#' @param horizon The forecast horizon, default to 1.
#' @param lags the number of consecutive lags in the model. 
#' @param fc.window Should the number of training observations be fixed or expanding when the sample increases?
#' @param fc.type Should the forecasts be recursive or direct?
#' @param ic Information criterion used to select the penalty parameter: BIC or AIC.
#' @param adaptive Optional, initial estimator for the adaptive Lasso. Can be chosen between: none (default, plain Lasso), ols, lasso, ridge. 
#' @param mc Parallelize the forecasts with the multicore package. If you want to parallelize the equation-by-equation esitmation of the lasso var, set mc=FALSE and mclas=TRUE. In both cases, use the optional parameter ncores to select the numbers of cores to use.
#' @param silent Should output be printed?
#' @param trend Should a linear trend be included in the model. 
#' @param post Optional, Should a post Lasso OLS be estimated, default FALSE.
#' @param ncores Optional, the number of cores to use in case of parallelization.
#' @param ... Optional options...
#'
#'
#' @return A list.
#' @export
forecast.lassovar<-function(dat,exo=NULL,fc.train,horizon=1,lags=1,fc.window=c('fix','expanding'),fc.type=c('recursive','direct'),ic=c('BIC','AIC'),adaptive=c('none','ols','lasso','group','ridge'),mc=FALSE,silent=FALSE,trend=FALSE,post=FALSE,ncores=NULL,... )
{
	# Matching a few args
	ic <- match.arg(ic)	
	fc.window <- match.arg(fc.window)
	fc.type <- match.arg(fc.type)
	adaptive <- match.arg(adaptive)
	
# Optional args
	argList	<- list(...)

# multicore ?
	if(!is.null(ncores))ncores<-ncores	else ncores	<-NULL

#Checking if exo exists and if so coercing to dataframe.
	if(!is.null(exo))exo<-as.data.frame(exo) else exo <-NULL	


# Storage init
	fc.lv	<-list('call'=match.call(),'err'=NULL,'coefficients'=NULL,'lambda'=NULL,'pred'=NULL)	

# Multicore?
	if(mc) if(is.null(ncores)){ncores=1}

# Number of forecasts
	nbr.fc	<-nrow(dat)-fc.train

# Some output 
if(!silent){
	cat('\n		-----------------------------		\n')	
	cat('Lassovar forecast\n')
	cat('Estimator: ',ifelse(is.null(argList$adaptive),'Lasso','Adaptive Lasso'),'\n',sep='')
	if(!is.null(argList$adaptive))cat('Initial Estimator: ',argList$adaptive,'\n',sep='')
	cat('Number of equations: ',ncol(dat),'\n',sep='')
	if(mc)cat('Number of cores used:',ncores)
	cat(fc.window,' window forecasts\n',sep='')
	cat(horizon,'-steps ahead ',fc.type,' forecasts. Initial training sample: ',fc.train,' observations.\n','Number of forecasts: ',nbr.fc,'\n',sep='')
	if(mc) cat('Forecast level multicore enabled, #cores: ',ncores,'\n',sep='')
	if(!is.null(argList$mclas))if(argList$mclas)cat('Equation level multicore enabled, #cores: ',ncores,'\n',sep='')
	cat('\n		-----------------------------		\n')	
}

# Progress bar here?

	if(!mc)fc.lassovar<-lapply(0:(nbr.fc-horizon),.fc.loop.lassovar
                             ,dat,fc.train,fc.window,lags
                             ,horizon,ic,exo=exo
							 ,adaptive=argList$adaptive
                             ,post=post,silent=silent
							 ,trend=trend)

	if(mc)fc.lassovar<-mclapply(0:(nbr.fc-horizon),.fc.loop.lassovar
                              ,dat,fc.train,fc.window,lags
                              ,horizon,ic,mc.cores=ncores,exo=exo
							  ,adaptive=argList$adaptive
                              ,post=post,silent=silent
							  ,trend=trend)
	

	fc.lv$coefficients<-list()
	if(post)fc.lv$post <- list()

	fc.lv$spectest <- array(NA,dim=c(3,ncol(dat),nbr.fc)
                          ,dimnames=list('Test'=c('Autocorr','Normality','Rsquared')
                                         ,'Equation'=1:ncol(dat),'Forecast'=1:nbr.fc))
	# Counter
	fc.i <- 1
	#aggregation
	for(lv.tmp in fc.lassovar)
	{
		fc.lv$lambda<-rbind(fc.lv$lambda,lv.tmp$lambda)
		fc.lv$err<-rbind(fc.lv$err,lv.tmp$fcerr)
		fc.lv$pred<-rbind(fc.lv$pred,lv.tmp$pred)
	 	fc.lv$coefficients <- c(fc.lv$coefficients,lv.tmp$coefficients)
		if(post) fc.lv$post<-c(fc.lv$post,lv.tmp$post)
	 	fc.lv$spectest[,,fc.i] <- lv.tmp$spectest
		fc.i <- fc.i+1
		
	}
return(fc.lv)
}


