\name{summary.lassovar}
\alias{summary.lassovar}

\title{
Returns the summary of a lassovar object.
}
\description{
	Prints a summary of a lassovar object.
}
\usage{
\method{summary}{lassovar}(object,...)
}

\arguments{
  \item{object}{An object of class lassovar.}
  \item{...}{An optional argument ic with value BIC or AIC. Default is BIC.}

}
\details{
	Prints a number of characteristics of the fitted model selected by criterion ic.
}
\value{nothing}
\author{
	Laurent Callot (lcallot@vu.nl)
}

\seealso{
	'lassovar' and 'predict.lassovar'.
}
\examples{
	data(LNprices)
	LNprices<-ts(LNprices,start=c(1964,3),frequency=12)
	las.mod	<-lassovar(LNprices,horizon=1,lags=1,criterion='BIC',type='equation',adaptive='ols')	
	summary(las.mod)

}
