# This file contains the two functions called by mclapply in forecast.lassovar.
# 
# 
# This function subsets the data into training sample and 'true' value of the forecast, calls the next.
.fc.loop.lassovar <-
function(fc,dat,fc.train,fc.window,lags,horizon,ic,exo,rest,post,adaptive,silent,trend,...)
{

	# begining of the training sample
	start.train	<- 1+ifelse(fc.window=='fix',fc,0)
	# end of the training sample.
	end.train	<- start.train + fc.train - 1 + ifelse(fc.window=='expanding',fc,0)
	
	# The training data
	train.dat <- dat[start.train:end.train,]
	if(!is.null(exo))train.exo <- as.matrix(exo[start.train:end.train,])
	else train.exo <-NULL
	
	# The 'true' forecast data 
	fc.dat		<-dat[end.train + 1,]	
	# Estimation
	fc.time		<-Sys.time()	
	fc.err		<-.fc.lassovar(train.dat,fc.dat,lags,horizon,ic,train.exo,rest,post,adaptive,trend,...)
	
	dif.time<-difftime(Sys.time(),fc.time,units='mins')
	if(!silent) cat('fc ',fc,' completed in ',round(dif.time,2),' minutes.\n',sep='')
	return(fc.err)
}

# This function estimates a lassovar, calls predict to get forecasts, computes the forecast errors and returns a list with forecast errors, predictions and diagnostics. 
.fc.lassovar<-function(train.dat,fc.dat,lags,horizon,ic='BIC',train.exo,rest,post,adaptive,trend,...)
{
	
	argList	<- list(...)
	mclas	<-argList$mclas
	if(is.null(mclas))mclas<-FALSE
	
	if(!is.null(argList$ncores))ncores <- ifelse(mclas,argList$ncores,1)
	else ncores <- 1
	
	y.var		<-.mkvar(train.dat,lags,horizon,train.exo,trend)
	
	las.mod		<-lassovar(dat=train.dat,horizon=horizon,lags=lags,ic=ic,exo=train.exo,mc=mclas,ncores=ncores,post=post,adaptive=adaptive,trend=trend)

	lv		<-list()
	
	lv$pred 	<-predict(las.mod,tail(y.var$x,1),ic=ic)
	# Forcast error
	lv$fcerr	<-(fc.dat-lv$pred)
	# Extracting and sparsifying the coefficients.
	lv$coefficients	<-Matrix(las.mod$coefficients,sparse=TRUE)
	if(!is.null(argList$post)){if(argList$post)lv$post<-Matrix(las.mod$post,sparse=TRUE)}
	# Saving the penalty parameters and the specification tests. 
	lv$lambda	<-las.mod$lambda
	lv$spectest	<-las.mod$spectest

	if(post) lv$post <- las.mod$post
	
	return(lv)
}


