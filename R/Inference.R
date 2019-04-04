
###################################
# General Loop
loop.iterative.procedure = function(Data,lyear,lmin=1,Kmax,Used.function,threshold=0.05){
  result  = list()
  n       = dim(Data)[1]
  SSwg    = c()
  loglik  = c()
  LogLg   = c()

  for (k in (1:Kmax)){
    res.segfunct=c()
    request=paste(paste0("res.segfunct=",Used.function,'(Data,k,lmin,lyear,threshold)'),sep="")
    eval(parse(text=request))
    SSwg[k]=res.segfunct$SSwg
    loglik[k]=-(n/2)*(log(2*pi))-(1/2)*(sum(log(res.segfunct$var.est.t)))-(1/2)*SSwg[k]
    seg       = matrix(res.segfunct$Tmu$end,ncol=k,nrow=1)
    LogLg[k] =  apply(seg,1,FUN=function(z) sum(log(diff(c(0,z))[diff(c(0,z))>0])))
  }
  result$SSwg=SSwg
  result$loglik=loglik
  result$LogLg=LogLg
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
