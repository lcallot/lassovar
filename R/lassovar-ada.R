#' @importFrom grpreg grpreg
#' @importFrom biglm biglm



# Group lasso adaptive weights
# Mostly Deprecated
.ada.grp.weights <-
function(y,x,ada,trend)
{
	nbr.lags	<-ncol(x)/ncol(y)
	grp.ind		<-sort(rep(1:nbr.lags,ncol(y)))
	grp.coef	<-NULL
	grp.coef.post	<-NULL

	greg.coef	<-NULL
	greg.coef.post	<-NULL

	gamma		<-1

	cat('Group Lasso Computation')
	cat('\n Equation: ',sep='')
	for(i in 1:ncol(y)){
		cat(i,' - ',sep='')

		greg.fit	<-grpreg(y=y[,i],X=x,group=grp.ind,family='gaussian',max.iter=500,eps=0.1)
		greg.bic	<-2*logLik(greg.fit)+log(nrow(y))*greg.fit$df
		lamb.min	<-greg.fit$lambda[which(greg.bic==min(greg.bic[-1]))]
		greg.pred	<-predict(greg.fit,X=x,lambda=lamb.min)
		greg.ipar	<-greg.fit$beta[,which(greg.bic==min(greg.bic[-1]))]
		greg.coef	<-cbind(greg.coef,greg.ipar[-1])

		greg.ipa.post	<-rep(0,length(greg.ipar)-1)
		greg.ipa.post[greg.ipar[-1]!=0]<-lm.fit(x=x[,greg.ipar[-1]!=0],y=y[,i])$coefficients
		greg.coef.post	<-cbind(greg.coef.post,greg.ipa.post)

	}
	rownames(greg.coef)=colnames(x)
	colnames(greg.coef)=colnames(y)
	rownames(greg.coef.post)=colnames(x)
	colnames(greg.coef.post)=colnames(y)
	cat('Group lasso ---- done','\n',sep='')

	grp.w<-list('ada'=ada,'b'=greg.coef,'b.post'=greg.coef.post,'w'=abs(greg.coef)^(-gamma),'gamma'=gamma)
	return(grp.w)
}




.ada.ols.weights <-
function(y,x,ada,mc,ncores)
{
	ada.w	<-list('ada'=ada)
	gamma	<-1

    if(!mc)	ols.b	<-lapply(1:ncol(y),function(i,y,x){lm(y[,i]~x)$coefficients},y,x)
    if(mc)	ols.b	<-mclapply(1:ncol(y),function(i,y,x){lm(y[,i]~x)$coefficients},y,x,mc.cores=ncores)
	
	ols.b	<- do.call(cbind,ols.b)

	ols.w	<-abs(ols.b[-1,])^(-gamma)
	ada.w$b	<-ols.b
	ada.w$w	<-ols.w
	ada.w$gamma	<-gamma

	return(ada.w)
}



# Lasso adaptive weights using the parameters from Lassovar
.ada.las.weights<-function(y,x,ada,ic='BIC',mc=NULL,ncores,dfmax,trend)
{
	ada.w	<-list('ada'=ada)
	gamma	<-1
	if(is.null(mc))mc<-FALSE

	lv.las	<-.lassovar.eq(y=y,x=x,degf.type=NULL,ada.w=NULL,ic=ic,mc=mc,ncores=ncores,alpha=1,dfmax=dfmax,trend=trend)
	
	ada.w$b	<-lv.las$coefficients
	ada.w$w	<-abs(lv.las$coefficients[-1,])^(-gamma)
	ada.w$gamma	<-gamma
	
	return(ada.w)
}




# Ridge regression for adaptive weights. 
.ada.ridge.weights<-function(y,x,ada,ic='BIC',mc=NULL,ncores,dfmax,trend)
{
	ada.w	<-list('ada'=ada)
	gamma	<-1
	if(is.null(mc))mc<-FALSE

	lv.las	<-.lassovar.eq(y=y,x=x,degf.type=NULL,ada.w=NULL,ic=ic,mc=mc,ncores=ncores,alpha=0,dfmax=dfmax,trend=trend)
	
	ada.w$b	<-lv.las$coefficients
	ada.w$w	<-abs(ada.w$b[-1,])^(-gamma)
	ada.w$gamma	<-gamma
	
	return(ada.w)
}

# A small function to compute the degrees of freedom for the ridge regresson
.ridge.df<-function(x,lambda){

	x<-(cbind(1,x))
	df.ridge <- NULL
	for(l in lambda){
		df.ridge<-c(df.ridge,sum(diag(x%*% solve(t(x) %*% x + l*diag(ncol(x)))%*% t(x))))
	}
	return(df.ridge)
}
