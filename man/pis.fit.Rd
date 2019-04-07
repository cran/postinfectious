\name{pis.fit}
\alias{pis.fit}
\title{Estimating the incubation period distribution of a post-infectious syndrome}
\description{
This function estimates the incubation period distribution of a post-infectious syndrome with maximum likelihood estimation. The incubation period distribution of the antecedent infection and the post-infectious syndrome are allowed to be lognormal ("\code{LN}"), Weibull ("\code{WB}") or gamma ("\code{GM}") distributed. The data set is allowed to have cases with the actecedent diseases whose incuation periods come from different distributions (see Examples).
}
\usage{
pis.fit(data,postinfect=c("LN","WB","GM"),theta)
}
\arguments{
 \item{data}{A data.frame containing at least 4 columns. The first two columns represent (1) the time between the symptom onset of the antecedent infection and post-infectious syndrome and (2) the incubation period distribution of the antecedent infection (only "\code{LN}", "\code{WB}" and "\code{GM}"). The last two columns refer to the parameters of the incubation period distribution of the antecedent infection; for "\code{LN}", they are \code{meanlog} and \code{sdlog} as in \code{dlnorm}; for "\code{WB}", they are \code{shape} and \code{scale} as in \code{dweibull}; for "\code{GM}", they are \code{shape} and \code{rate} as in \code{dgamma}.}
 \item{postinfect}{The incubation period distribution of the post-infectious disease. It can only be "\code{LN}", "\code{WB}" and "\code{GM}".}
 \item{theta}{A vector of two numbers as the initial value for optimisation.}
}
\details{
For each observed case, let \eqn{S_{0}} and \eqn{S} be the incubation period of the antecedent infection and post-infectious syndrome, respectively. As the antecedent infection is the antigenic factor of the post-infectious syndrome, they both share the same time of infection exposure. The difference between \eqn{S_{0}} and \eqn{S}, denoted by \eqn{X}, is the time between the two symptom onsets. Also let \eqn{\theta_{0}} and \eqn{\theta} be the set of the parameters of the distribution of \eqn{S_{0}} and \eqn{S} then the likelihood of such observed case is given by,
\deqn{\int_{-\infty}^{\infty}f_0(S_0,\theta_0)f(S_0+X,\theta)dS_0}
where \eqn{f_0} and \eqn{f} are the probability density function of \eqn{S_{0}} and \eqn{S}, respectively. \eqn{\theta} is then estimated by maximising the sum of likelihood of all observed cases.
}
\value{
  \item{Parameter}{Estimates of the parameters of the incubation period distribution of the post-infectious syndrome.}
  \item{SE}{Standard errors of \code{Parameter}}
  \item{AIC}{Akaike Information Criterion.}
  \item{Convergence}{The convergence message of \code{optim}}
  \item{Median}{The median incubation period distribution of the post-infectious syndrome.}
  \item{Theta.initial}{Initial values used in \code{optim}}
  \item{Distribution}{The Distribution assumed in the estimation, i.e. "\code{LN}", "\code{WB}" or "\code{GM}".}
}
\author{Char Leung}
\examples{
#generate artificial data
S<-c(56,37,32,7,8,3,5)
S0<-c(2,1,3,1,1,1,3)
X<-S-S0
f0<-c(rep("LN",4),rep("WB",3))
phi<-matrix(c(rep(c(0,1),4),rep(c(1,2),3)),byrow=TRUE,ncol=2)
data<-data.frame(X,f0,phi)
pis.fit(data,"LN",theta=c(2.5,1))
}

