\name{pis.fit.boots}
\alias{pis.fit.boots}
\title{Bootstrap estimates of the output in \code{pis.fit}}
\description{
This function creates bootstrap estimates of the output of \code{pis.fit} by creating bootstrap samples
}
\usage{
pis.fit.boots(data,postinfect=c("LN","WB","GM"),theta,n.boots=1000,collective=100)
}
\arguments{
\item{data}{A data.frame containing at least 4 columns. The first two columns represent (1) the time between the symptom onset of the antecedent infection and post-infectious syndrome and (2) the incubation period distribution of the antecedent infection (only "\code{LN}", "\code{WB}" and "\code{GM}"). The last two columns refer to the parameters of the incubation period distribution of the antecedent infection; for "\code{LN}", they are \code{meanlog} and \code{sdlog} as in \code{dlnorm}; for "\code{WB}", they are \code{shape} and \code{scale} as in \code{dweibull}; for "\code{GM}", they are \code{shape} and \code{rate} as in \code{dgamma}.}
\item{postinfect}{The incubation period distribution of the post-infectious disease. It can only be "\code{LN}", "\code{WB}" and "\code{GM}".}
\item{theta}{Text input only and it is an R expression to be evaluated (i.e. \code{eval}) so as to create initial values used in \code{optim}. The reason behind this is to allow random numbers as the initial values in optimisation. See Examples.}
\item{n.boots}{The number of bootstrap samples.}
\item{collective}{The number of bootstrap samples to be estimated at once as the estimation process uses the \code{apply} function.}
 }
\value{
Same as those in \code{pis.fit}.
}
\author{Char Leung}
\examples{
S<-c(56,37,32,7,8,3,5)
S0<-c(2,1,3,1,1,1,3)
X<-S-S0
f0<-c(rep("LN",4),rep("WB",3))
phi<-matrix(c(rep(c(0,1),4),rep(c(1,2),3)),byrow=TRUE,ncol=2)
data<-data.frame(X,f0,phi)
pis.fit.boots(data,"LN",theta="c(runif(1,2,3),runif(1,0,1))",n.boots=20,collective=15)
}
\seealso{pis.fit}
