
###################################
######## Inference procedure for a fixed K
Seg_funct_totK <-function(Data,var.est.month,K,lmin,lyear,threshold,tol){
  var.est.t=var.est.month[as.numeric(Data$month)]
  
  period=periodic_estimation_tot_init(Data,var.est.t,lyear)
  auxiliar_data <- Data
  auxiliar_data$signal=Data$signal-period$predict
  segmentation=SegMonthlyVarianceK(auxiliar_data,K,lmin,var.est.t)

  maxIter = 100
  Diff    = 2*tol
  Iter=0

  while ((Diff  > tol) & (Iter < maxIter))
  {
    Iter = Iter +1
    auxiliar_data$signal=Data$signal-segmentation$mean.est.t
    periodi=periodic_estimation_tot(auxiliar_data,var.est.t,lyear)

    auxiliar_data$signal=Data$signal-periodi$predict
    segmi=SegMonthlyVarianceK(auxiliar_data,K,lmin,var.est.t)

    if (Iter == 2)
    {
      t2 = c(period$predict,segmentation$mean.est.t)
    }
    if (Iter == 3)
    {
      t1 = c(period$predict,segmentation$mean.est.t)
      t0 = c(periodi$predict,segmi$mean.est.t)
      tp0 = (t2-t1)/sum((t1==t2)+((t2-t1)^2)) + (t0-t1)/sum((t1==t0)+((t0-t1)^2))
      tp0 = t1 + tp0 / sum((tp0==0) + tp0^2)
    }
    if (Iter > 3)
    {
      t2 = t1
      t1 = t0
      t0 = c(periodi$predict,segmi$mean.est.t)
      tp1 = tp0
      tp0 = (t2-t1)/sum((t1==t2)+((t2-t1)^2)) + (t0-t1)/sum((t1==t0)+((t0-t1)^2))
      tp0 = t1 + tp0 / sum((tp0==0) + tp0^2)
      Diff = sum((tp0-tp1)^2)
    }

    period=periodi
    segmentation=segmi
  }

  
  segmK=c()
  segmK$Tmu   = segmentation$Tmu
  segmK$SSwg  = segmentation$SSwg[K]
  segmK$LogLg = segmentation$LogLg[K]
  segmK$f     = period$predict
  segmK$coeff = period$coeff
  return(segmK)

}


Seg_funct_selbK <-function(Data,var.est.month,K,lmin=1,lyear,threshold,tol){
  var.est.t=var.est.month[as.numeric(Data$month)]
  
  period=periodic_estimation_selb_init(Data,var.est.t,lyear,threshold)
  auxiliar_data <- Data
  auxiliar_data$signal=Data$signal-period$predict
  segmentation=SegMonthlyVarianceK(auxiliar_data,K,lmin,var.est.t)

  maxIter = 100
  Diff  = 2*tol
  Iter  = 0

  while ((Diff  > tol) & (Iter < maxIter))
  {
    Iter = Iter +1
    auxiliar_data$signal=Data$signal-segmentation$mean.est.t
    periodi=periodic_estimation_selb(auxiliar_data,var.est.t,lyear,threshold)

    auxiliar_data$signal=Data$signal-periodi$predict
    segmi=SegMonthlyVarianceK(auxiliar_data,K,lmin,var.est.t)

    if (Iter == 2)
    {
      t2 = c(period$predict,segmentation$mean.est.t)
    }
    if (Iter == 3)
    {
      t1 = c(period$predict,segmentation$mean.est.t)
      t0 = c(periodi$predict,segmi$mean.est.t)
      tp0 = (t2-t1)/sum((t1==t2)+((t2-t1)^2)) + (t0-t1)/sum((t1==t0)+((t0-t1)^2))
      tp0 = t1 + tp0 / sum((tp0==0) + tp0^2)
    }
    if (Iter > 3)
    {
      t2 = t1
      t1 = t0
      t0 = c(periodi$predict,segmi$mean.est.t)
      tp1 = tp0
      tp0 = (t2-t1)/sum((t1==t2)+((t2-t1)^2)) + (t0-t1)/sum((t1==t0)+((t0-t1)^2))
      tp0 = t1 + tp0 / sum((tp0==0) + tp0^2)
      Diff = sum((tp0-tp1)^2)
    }

    period=periodi
    segmentation=segmi

  }

  
  segmK=c()
  segmK$Tmu   = segmentation$Tmu
  segmK$SSwg  = segmentation$SSwg[K]
  segmK$LogLg = segmentation$LogLg[K]
  segmK$f     = period$predict
  segmK$coeff = period$coeff
  return(segmK)

}


###################################
######## Robust estimation of the variances
RobEstiMonthlyVariance <- function(Y){
  Kmonth <- length(unique(Y$month))
  z <- stats::aggregate(signal ~ month + year, data = Y, diff)
  sigma.est <- sapply(1:Kmonth,function(i) {
    e <- subset(z,z$month==levels(z$month)[i])
    ee <- unlist(e$signal)
    robustbase::Qn(ee, constant = 1 / (sqrt(2) * stats::qnorm(5/8)))/sqrt(2)

  })
  return(sigma.est)
}




###################################
######## Functions for segmentation
SegMonthlyVarianceK=function(Data,K,lmin,var.est.t){
  result=list()
  vh=3
  matD  = Gsegmentation(Data,vh,lmin,var.est.t)
  out   = DynProg(matD,K)

  Tmu=c()
  Tmu=FormatOptSegK(out$t.est[K,1:K],Data,var.est.t)
  mean.est.t  = rep(Tmu$mean,diff(c(0,Tmu$end)))

  result$Tmu= Tmu
  result$res.LoopK=out
  result$mean.est.t = mean.est.t
  result$SSwg=out$J.est
  result$LogLg=apply(out$t.est,1,FUN=function(z) sum(log(diff(c(0,z))[diff(c(0,z))>0])))

  return(result)
}

# Cost Matrix
Gsegmentation<-function(x,vh,lmin,var.t=NULL) {
  # segmentation with homogeneous variance
  if (vh==1){
    x=x$signal
    n = length(x)
    matD=matrix(Inf,n,n)
    x2=x^2
    x2i=cumsum(x2)
    xi=cumsum(x)
    x2i=x2i[lmin:n]
    xi=xi[lmin:n]
    matD[1,lmin:n]=x2i-((xi^2)/(lmin:n))
    nl=n-lmin+1

    for (i in 2:nl){
      ni=n-i-lmin+3
      x2i=x2i[2:ni]-x2[i-1]
      xi=xi[2:ni]-x[i-1]
      deno=((i+lmin-1):n)-i+1
      matD[i,(i+lmin-1):n]=x2i-((xi^2)/deno)
    }
  }

  # segmentation with heterogeneous variances
  ## by segment
  if (vh==2) {
    x=x$signal
    n = length(x)
    matD=matrix(Inf,n,n)
    x2=x^2
    x2i=cumsum(x2)
    xi=cumsum(x)
    x2i=x2i[lmin:n]
    xi=xi[lmin:n]
    matD[1,lmin:n]=(lmin:n)*log(( x2i-(xi^2)/(lmin:n) ) / (lmin:n) )
    nl=n-lmin+1

    for (i in 2:nl){
      ni=n-i-lmin+3
      x2i=x2i[2:ni]-x2[i-1]
      xi=xi[2:ni]-x[i-1]
      matD[i,(i+lmin-1):n]=(lmin:(n-i+1))*(log((x2i-(xi^2)/(lmin:(n-i+1)))/(lmin:(n-i+1))))
    }
  }


  ## time-dependant
  if (vh==3) {

    x=x$signal
    b=var.t
    n = length(x)
    #print(n)
    matD=matrix(Inf,n,n)
    x2=x^2/b
    x2i=cumsum(x2)
    xb=x/b
    xbi=cumsum(xb)
    sig=1/b
    sigi=cumsum(sig)
    matD[1,lmin:n]=x2i-(xbi^2)/sigi
    nl=n-lmin+1

    for (i in 2:nl) {
      ni=n-i-lmin+3
      x2i=x2i[2:ni]-x2[i-1]
      xbi=xbi[2:ni]-xb[i-1]
      sigi=sigi[2:ni]-sig[i-1]
      matD[i,(i+lmin-1):n]=x2i-(xbi^2)/sigi
    }
  }
  invisible(matD)
}

# Dynamic Programming
DynProg<-function(matD,Kmax)
{
  N<-dim(matD)[1]
  if (Kmax>N){cat("Kmax ", Kmax, "is greater than N ",N,"\n")
    cat("Kmax is supposed to be equal to N :", N,"\n")
    Kmax <- N/2
  }

  I<-matrix(Inf,Kmax,N)
  t<-matrix(0,Kmax,N)
  I[1,]=matD[1,]
  matD=t(matD)


  if (Kmax>2)
  {
    for (k in 2:(Kmax-1)){
      for (L in k:N)
      {
        I[k,L]<-min(I[(k-1),1:(L-1)]+matD[L,2:L])
        if(I[k,L]!=Inf){
          t[k-1,L]<-which.min(I[(k-1),1:L-1]+matD[L,2:L])
        }
      }
    }
    I[Kmax,N]<-min(I[Kmax-1,1:(N-1)]+matD[N,2:N])
    if(I[Kmax,N]!=Inf){
      t[Kmax-1,N]<-which.min(I[(Kmax-1),1:N-1]+matD[N,2:N])
    }
  } else if (Kmax==2) {

    I[Kmax,N]<-min(I[Kmax-1,1:(N-1)]+matD[N,2:N])
    if(I[Kmax,N]!=Inf){
      t[Kmax-1,N]<-which.min(I[(Kmax-1),1:N-1]+matD[N,2:N])
    }
  }
  # *** Calcul des instants de ruptures ***

  t.est<-matrix(0,Kmax,Kmax)
  diag(t.est)<-N

  if (Kmax>=2) {

    for (K in 2:Kmax){

      for (k in seq(K-1,1,by=-1))
      {
        if(t.est[K,k+1]!=0){
          t.est[K,k]<-t[k,t.est[K,k+1]]
        }
      }
    }
  }
  list(J.est = I[,N],t.est = t.est)
}

# Format of the segmentation result
FormatOptSegK <- function(breakpointsK,Data,v){
  K     = length(breakpointsK)
  rupt  = matrix(Inf,ncol = 2 , nrow= K)
  bp    = breakpointsK
  rupt[,2]  = bp
  bp        = bp +1
  rupt[,1]  = c(1, bp[1:K-1])

  Data.var   = Data$signal/v
  mean.est.k = apply(rupt,1,FUN=function(z) sum(Data.var[z[1]:z[2]],na.rm=TRUE))
  var.est.k = apply(rupt,1,FUN=function(z) sum(1/(v[z[1]:z[2]]),na.rm=TRUE))
  mean.est.k=mean.est.k/var.est.k
  Tmu=data.frame(rupt,mean.est.k)
  colnames(Tmu) = c("begin","end","mean")
  return(Tmu)
}

###################################
######## Functions for functional
periodic_estimation_tot=function(Data,var.est.t,lyear){
  DataF=Data
  DataF$t=c(as.numeric(DataF$date-DataF$date[1]))/86400
  for (i in 1:4){
    cosX=cos(i*DataF$t*(2*pi)/lyear)
    sinX=sin(i*DataF$t*(2*pi)/lyear)
    DataF=cbind(DataF,cosX,sinX)
    colnames(DataF)[(dim(DataF)[2]-1):(dim(DataF)[2])]=c(paste0('cos',i),paste0('sin',i))
  }
  reg=stats::lm(signal~-1+cos1+sin1+cos2+sin2+cos3+sin3+cos4+sin4,weights=1/var.est.t,data=DataF)
  coeff=base::summary(reg)$coefficients[,1]
  result=list()
  result$predict=stats::predict(reg,DataF)
  result$coeff=coeff
  return(result)
}


periodic_estimation_tot_init=function(Data,var.est.t,lyear){
  DataF=Data
  DataF$t=c(as.numeric(DataF$date-DataF$date[1]))/86400
  for (i in 1:4){
    cosX=cos(i*DataF$t*(2*pi)/lyear)
    sinX=sin(i*DataF$t*(2*pi)/lyear)
    DataF=cbind(DataF,cosX,sinX)
    colnames(DataF)[(dim(DataF)[2]-1):(dim(DataF)[2])]=c(paste0('cos',i),paste0('sin',i))
  }
  reg=stats::lm(signal~-1+cos1+sin1+cos2+sin2+cos3+sin3+cos4+sin4,data=DataF)
  coeff=base::summary(reg)$coefficients[,1]
  result=list()
  result$predict=stats::predict(reg,DataF)
  result$coeff=coeff
  return(result)
}


periodic_estimation_selb=function(Data,var.est.t,lyear,threshold=0.05){
  DataF=Data
  DataF$t=c(as.numeric(DataF$date-DataF$date[1]))/86400
  num.col=dim(DataF)[2]
  for (i in 1:4){
    cosX=cos(i*DataF$t*(2*pi)/lyear)
    sinX=sin(i*DataF$t*(2*pi)/lyear)
    DataF=cbind(DataF,cosX,sinX)
    colnames(DataF)[(dim(DataF)[2]-1):(dim(DataF)[2])]=c(paste0('cos',i),paste0('sin',i))
  }
  reg=stats::lm(signal~-1+cos1+sin1+cos2+sin2+cos3+sin3+cos4+sin4,weights=1/var.est.t,data=DataF)
  res.coeff=summary(reg)$coefficients
  names.coeff=rownames(res.coeff)
  #Selection of significant coefficients
  rg=which(res.coeff[,4]<threshold)

  if (length(rg)>=1){
    names.Selected=names.coeff[rg]
    n.Selected=length(names.Selected)
    DataFF=DataF[,c(1:num.col,which(colnames(DataF) %in%  names.Selected ))] 
    if (n.Selected >=2){
      a=names.Selected[1]
      for (i in 1:(n.Selected-1)){
      a=paste(a,names.Selected[i+1],sep="+")
      }
      } else {
      a=names.Selected[1]
      }
    request=paste(paste0("reg.Selected=stats::lm(signal~-1+",a,",weights=1/var.est.t,data=DataFF)"),sep="")
    eval(parse(text=request))
    pred=reg.Selected$fitted.values
    coeff=reg.Selected$coefficients
  } else {
    pred=rep(0,length(DataF$signal))
    coeff=0
    names(coeff)="no selected coeff"
  }
  result=list()
  result$predict=pred
  result$coeff=coeff
  return(result)
}
  
  

periodic_estimation_selb_init=function(Data,var.est.t,lyear,threshold=0.05){
  DataF=Data
  DataF$t=c(as.numeric(DataF$date-DataF$date[1]))/86400
  num.col=dim(DataF)[2]
  for (i in 1:4){
    cosX=cos(i*DataF$t*(2*pi)/lyear)
    sinX=sin(i*DataF$t*(2*pi)/lyear)
    DataF=cbind(DataF,cosX,sinX)
    colnames(DataF)[(dim(DataF)[2]-1):(dim(DataF)[2])]=c(paste0('cos',i),paste0('sin',i))
  }
  
  
  reg=stats::lm(signal~-1+cos1+sin1+cos2+sin2+cos3+sin3+cos4+sin4,data=DataF)
  res.coeff=summary(reg)$coefficients
  names.coeff=rownames(res.coeff)
  #Selection of significant coefficients
  rg=which(res.coeff[,4]<threshold)
  
  if (length(rg)>=1){
    names.Selected=names.coeff[rg]
    n.Selected=length(names.Selected)
    DataFF=DataF[,c(1:num.col,which(colnames(DataF) %in%  names.Selected ))] 
    if (n.Selected >=2){
      a=names.Selected[1]
      for (i in 1:(n.Selected-1)){
        a=paste(a,names.Selected[i+1],sep="+")
      }
    } else {
      a=names.Selected[1]
    }
    request=paste(paste0("reg.Selected=stats::lm(signal~-1+",a,",data=DataFF)"),sep="")
    eval(parse(text=request))
    pred=reg.Selected$fitted.values
    coeff=reg.Selected$coefficients
  } else {
    pred=rep(0,length(DataF$signal))
    coeff=0
    names(coeff)="no selected coeff"
  }
  result=list()
  result$predict=pred
  result$coeff=coeff
  return(result)
}
