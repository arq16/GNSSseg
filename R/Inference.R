
###################################
# General Loop
Loop.K.procedure = function(Data,var.est.month,lyear,lmin=1,Kmax,Used.function,threshold,tol){
  result  = list()
  n       = dim(Data)[1]
  #Iterative procedure
  result=lapply(1:Kmax, function(i){
    res.segfunct=c()
    request=paste(paste0("res.segfunct=",Used.function,'(Data,var.est.month,i,lmin,lyear,threshold,tol)'),sep="")
    eval(parse(text=request))
  })
  return(result)
}



###################################
######## Model selection criteria
MLcriterion<-function(J,Kseq,S=0.75)
{
  Kmax=length(Kseq)
  Jtild=(J[Kmax]-J)/(J[Kmax]-J[1])*(Kseq[Kmax]-Kseq[1])+1
  D=diff(diff(Jtild))

  if (length(which(D>=S))>0){
    Kh <- max(which(D>=S))+1
  }else {
    Kh <- 1
  }
  Kh=Kseq[Kh]
  return(Kh)
}

BMcriterion<-function(J,pen)
{
  Kmax=length(J)
  Kseq=1:Kmax
  k=1
  kv=c()
  dv=c()
  pv=c()
  dmax=1
  while (k<Kmax) {
    pk=(J[(k+1):Kmax]-J[k])/(pen[k]-pen[(k+1):Kmax])
    pm=max(pk)
    dm=which.max(pk)
    dv=c(dv,dm)
    kv=c(kv,k)
    pv=c(pv,pm)
    if (dm>dmax){
      dmax=dm
      kmax=k
      pmax=pm
    } #end
    k=k+dm
  } #end

  pv=c(pv,0)
  kv=c(kv,Kmax)
  dv=diff(kv);
  dmax=max(dv)
  rt=which.max(dv)
  pmax=pv[rt[length(rt)]]
  alpha=2*pmax
  km=kv[alpha>=pv]
  km=km[1]
  Kh =Kseq[km]
  return(Kh)
}

mBICcriterion <- function(SSwg,LogLg,n,Kseq){
  prova=list()
  Kh=c()
  mBIC=c()
  mBIC=(1/2)*(-SSwg)-(1/2)*LogLg+(3/2-Kseq)*log(n)
  Kh.mBIC=which.max(mBIC)
  prova$Kh=Kh.mBIC
  prova$mBIC=mBIC
  return(prova)
}
