#' @title Estimating the Incubation Period Distribution of Post-Infectious Syndrome
#'
#' @description This package estimates the incubation period distribution of post-infectious syndrome which is defined as the time between the symptom onset of the antecedent infection and that of the post-infectious syndrome.
#'
#' @param data,postinfect,theta,dataM,postinfect0,theta0,n.boots,collective
#'
#' @return sumL.i,sumL.i0,res
#'
#' @examples  pis.fit(med.data,"LN",c(runif(1,0,5),runif(1,0,1)))), pis.fit.boots(med.data,"LN",theta="c(runif(1,2,3),runif(1,0,1))",n.boots=100,collective=20)
#'
#' @export   pis.fit.boots,pis.fit


pis.fit<-function(data,postinfect=c("LN","WB","GM"),theta){

  if(sum(!is.numeric(data[,1]))>0) stop("Non-numeric in X")
  if(sum(!is.numeric(unlist(c(data[,3:4]))))>0) stop("Non-numeric in Theta0")
  if(sum(!(data[,2]%in%c("LN","WB","GM")))>0) stop("Undefined distributions in f0: can only be 'LN', 'WB' or 'GM'")

  dist.read<-matrix(c("LN","WB","GM","dlnorm","dweibull","dgamma","qlnorm","qweibull","qgamma"),nrow=3,ncol=3)
  f0.list<-data[,2:4];f0.list<-f0.list[!duplicated(f0.list),]
  L.j<-NULL
  L.i<-function(B,x0,antecedent.dist,postinfect.dist,theta0){
    eval(parse(text=paste("L.j<-function(z){-log(integrate(function(S0){",antecedent.dist,"(S0,",theta0[1],",",theta0[2],")*",postinfect.dist,"(z+S0,B[1],B[2])","},-Inf,Inf)$value)}",sep="")))
    return(sum(apply(t(t(x0)),1,L.j)))}

  L<-function(B,L.f0.list=f0.list,L.dist.read=dist.read,L.data=data){
    sumL.i<-0
    for(i in 1:nrow(L.f0.list)){
      x.i<-L.data[L.data[,2]==L.f0.list[i,1]&L.data[,3]==L.f0.list[i,2]&L.data[,4]==L.f0.list[i,3],1]
      sumL.i<-sumL.i+L.i(B,x0=x.i,antecedent.dist=L.dist.read[L.dist.read[,1]==L.f0.list[i,1],2],postinfect.dist=L.dist.read[L.dist.read[,1]==postinfect,2],theta0=L.f0.list[i,2:3])
    }
    return(sumL.i)
  }

  res<-optim(theta,L,hessian=TRUE,method="L-BFGS-B",control=list(maxit=30000))
  aic<-4+2*res$value
  vcov<-solve(res$hessian)
  se<-c(sqrt(diag(vcov)))
  para<-res$par
  convg<-res$convergence
  quan<-eval(parse(text=paste(dist.read[dist.read[,1]==postinfect,3],"(0.5,",para[1],",",para[2],")",sep="")))
  return(list(Parameter=para,SE=se,AIC=aic,Convergence=convg,Median=quan,Theta.initial=theta,Distribution=postinfect))}



pis.fit0<-function(dataM,postinfect0=c("LN","WB","GM"),theta0){
  if(!is.matrix(dataM)) stop("please supply matrix form")
  res0<-list(Parameter=rep(NA,2),SE=rep(NA,2),Convergence=NA,Median=NA,Theta.initial=rep(NA,2),AIC=NA)
  X<-dataM[,1];dist.label<-numeric(nrow(dataM));thetaM<-matrix(dataM[,3:4],ncol=2)
  for(i in 1:nrow(dataM)){
    if(dataM[i,2]==1){dist.label[i]<-"LN"}
    if(dataM[i,2]==2){dist.label[i]<-"WB"}
    if(dataM[i,2]==3){dist.label[i]<-"GM"}
  }
  data0<-data.frame(X,dist.label,thetaM)
  tryCatch(res0<-pis.fit(data=data0,postinfect=postinfect0,theta=theta0),error=function(e) NULL)
  sumL.i0<-unlist(list(Parameter=res0$Parameter,SE=res0$SE,AIC=res0$AIC,Convergence=res0$Convergence,Median=res0$Median,Theta.initial=res0$Theta.initial))
  return(sumL.i0)
}



pis.fit.boots<-function(data,postinfect=c("LN","WB","GM"),theta,n.boots=1000,collective=100){
  if(sum(!is.numeric(data[,1]))>0) stop("Non-numeric in X")
  if(sum(!is.numeric(unlist(c(data[,3:4]))))>0) stop("Non-numeric in Theta0")
  if(sum(!(data[,2]%in%c("LN","WB","GM")))>0) stop("Undefined distributions in f0: can only be 'LN', 'WB' or 'GM'")
  theta.eval<-eval(parse(text=theta))
  if(sum(!is.numeric(theta.eval))>0) stop("Undefined theta")
  q0<-0.5
  dataM<-data.matrix(data, rownames.force = NA)
  for(i in 1:nrow(dataM)){
    if(data[i,2]=="LN"){dataM[i,2]<-1}
    if(data[i,2]=="WB"){dataM[i,2]<-2}
    if(data[i,2]=="GM"){dataM[i,2]<-3}
  }
  N<-nrow(data)
  res<-NULL
  cat(sprintf("Simulating %s \n",""))
  pb <- txtProgressBar(min = 0, max = n.boots, style = 3,char="|")
  repeat{
    compiled.data<-array(NA,dim=c(N,4,collective))
    for(j in 1:collective){
      compiled.data[,,j]<-dataM[sample(1:N,N,replace=TRUE),]
    }
    theta.eval<-eval(parse(text=theta))
    res0<-suppressWarnings(apply(compiled.data,3,pis.fit0,postinfect0=postinfect,theta0=theta.eval))
    res0<-res0[,colSums(is.na(res0))<1 & res0[6,]==0]
    res<-cbind(res,res0)
    if(ncol(res)>n.boots){res<-res[,1:n.boots]}
    i<-ncol(res)
    setTxtProgressBar(pb, i)
    if(i>(n.boots-1)) break}
  close(pb)
  colnames(res)<-NULL
  return(res)
}

