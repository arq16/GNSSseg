#' Homogeneization of GNSS series
#'
#' fit a segmentation in the mean model by taken into account for a functional part (estimated with a Fourier decomposition) and a heterogeneous variance (monthly variance)
#'
#' @param Data a data frame, with size [n x 2], containing the signal (e.g. the daily GPS-ERAI series for GNSS) and the dates (in format yyyy-mm-dd of type "calendar time" (class POSIXct))
#' @param lyear the length of the year in the signal. Defalut is 365.25
#' @param lmin the minimum length of the segments. Default is 1
#' @param Kmax the maximal number of segments (must be lower than n). Default is 40#[n/2 log(n)]
#' @param selection.K a name indicating the model selection criterion to select the number of segments K (\code{mBIC}, \code{Lav}, \code{BM_BJ} or \code{BM_slope}). \code{"none"} indicates that no selection is claimed and the procedure considers \code{Kmax} segments or \code{Kmax}-1 changes. If \code{selection.K="All"}, the results for the four possible criteria are given. Default is \code{"mBIC"}
#' @param S the threshold used in the Lav criterion. Default is 0.75
#' @param f a boolean indicating if the functional part is taking into account for in the model. Default is FALSE
#' @param selection.f a boolean indicating if a selection on the functions of the Fourier decomposition of order 4 is performed. Default is FALSE
#' @param threshold a numeric value lower than 1. Default is 0.05
#' @param tol the stopping rule for the iterative procedure. Default is 1e-4
#'
#' @return A file containing
#' \itemize{
#' \item \code{selected.K} that corresponds to the selected number of segments. If \code{selection.K="none"}, the number of segments is fixed to \code{Kmax}
#' \item \code{segmentation} that corresponds to the estimation of the segmentation parameters (the begin and the end positions of each segment with the estimated mean)
#' \item \code{functional} that corresponds to the estimation of the functional part
#' \item \code{loglik} that corresponds to the log-likelihood for k=1,...,\code{Kmax}. If \code{selection.K="none"}, it contains only the log-likelihood for \code{Kmax} segments
#' \item \code{variances} that corresponds to the estimated variances for each fixed interval (e.g. the months).
#' \item \code{mBIC} that corresponds to the values of the mBIC criterion for k=1,...,\code{Kmax} if it is required (\code{selection.K="mBIC"} or \code{selection.K="All"})
#' }
#' If \code{selection.K="All"}, the outputs \code{selected.K}, \code{segmentation} and \code{functional} are each a list containing the corresponding result for the four model selection criteria
#'
#' @details
#' The function performs homogeneization of GNSS series. The considered model is a segmentation in the mean model (detection of abrupt changes model) in which a functional part is added and with heterogeneous variances on fixed intervals. On GNSS series, the latter intervals corresponds to the month.
#' The inference procedure consists for a fixed number of segments K in first estimating robustly the variances and second estimating the segmentation and the functional parts using an iterative procedure. The solution is obtained for k=1,...,\code{Kmax} segments. Then the "best k" is chosen using model selection criteria. The possible criteria are \code{mBIC} the modified BIC criterion REFEREF, \code{Lav} the criterion proposed by REFEF, \code{BM_BJ} and \code{BM_slope} the criteria proposed by REFEF where the penalty constant is calibrated using the Biggest Jump and the slope respectively REFERF.
#' \itemize{
#' \item The data is a data frame with 2 columns: $signal is the signal to be homogeneized (a daily series) and $date is the date. The date will be in format yyyy-mm-dd of type "calendar time" (class POSIXct).
#' \item The function part is estimated using a Fourier decomposition of order 4 with \code{selection.f=FALSE}. \code{selection.f=TRUE} consists in selecting twice functions of the Fourier decomposition of order 4, first by using a stepwise procedure using AIC and second by considering among them the signficative ones (for which p.values are lower than \code{threshold})
#' \item If \code{selection.K="none"} the procedure is performed with \code{Kmax} segments.
#' \item Missing data in the signal are accepted.
#' }
#' Note that if \code{f=FALSE}, only a segmentation is performed.
#'
#' @examples
#' data(Data)
#' lyear=365.25
#' Kmax=8
#' lmin=1
#' result=GNSSseg(Data,lyear,Kmax=Kmax)
#' plot_GNSS(Data,result$segmentation,result$functional)
#' @export

GNSSseg=function(Data,lyear=365.25,lmin=1,Kmax=40,selection.K="mBIC",S=0.75,f=TRUE,selection.f=FALSE,threshold=0.05,tol=1e-4){
  result  = list()
  Data.X  = c()
  cond1=TRUE
  cond2=TRUE
  n.present = length(which(!is.na(Data$signal)))

  #The conditions to be fulfilled
  if (class(Data$date)[1]!="POSIXct"){
    cond1=FALSE
    cat("date must be in format yyyy-mm-dd of type GMT in class POSIXct/POSIXt")
    }
  if (Kmax >n.present) {
    cond2=FALSE
    cat("The maximal number of segments Kmax", Kmax," needs to be lower than the length of the series without NA that is " ,n.present,"\n")}
  if ((cond1==TRUE) & (cond2==TRUE)){
      Data$year=format(Data$date,format='%Y')
      Data$month=format(Data$date,format='%m')
      Data$year=as.factor(Data$year)
      Data$month=as.factor(Data$month)


      #For NA
      present.data = which(!is.na(Data$signal))
      missing.data = which(is.na(Data$signal))
      Data.X       = Data[present.data,]
      Data.X$month = droplevels(Data.X$month)
      Data.X$year  = droplevels(Data.X$year)

      #Kmax
      n.Data=length(Data[,1])
      n.X=length(Data.X[,1])
      if (is.null(Kmax)){
        Kmax=floor(n.X/(2*log(n.X)))
      }
      Kseq=1:Kmax


      #Used function for NA
      add_NA=function(res,n.Data,present.data){
        res.with.NA=list()
        seg = res$Tmu$end
        f   = res$f
        Tmu = res$Tmu

        seg.Data=seg
        seg.Data=present.data[seg]

        f.Data=rep(NA,n.Data)
        f.Data[present.data]=f

        if (length(seg.Data)==1){
          Tmu$end=n.Data
        } else {
          Tmu$end=seg.Data
          Tmu$begin=c(0,seg.Data[1:(length(seg.Data)-1)])+1
        }
        res.with.NA$Tmu=Tmu
        res.with.NA$f=f.Data
        return(res.with.NA)
      }
      #

      if (f==TRUE){
        #Option for estimating f
        Used.function=c()
        if (selection.f==TRUE){
          Used.function='Seg_funct_selbK'
          } else{
            Used.function='Seg_funct_totK'
            }

        if (selection.K=="none"){
          res.segfunct=c()
          request=paste(paste0("res.segfunct=",Used.function,'(Data.X,Kmax,lmin,lyear,threshold,tol)'),sep="")
          eval(parse(text=request))
          Tmu=c()
          fh=c()
          varh=c()
          res.segfunct.with.NA<- add_NA(res.segfunct,n.Data,present.data)
          Tmu=res.segfunct.with.NA$Tmu
          fh=res.segfunct.with.NA$f
          varh=res.segfunct$sigma.est.month^2
          Kh=Kmax
          loglik=-(n.X/2)*(log(2*pi))-(1/2)*(sum(log(res.segfunct$var.est.t)))-(1/2)*res.segfunct$SSwg
          mBIC=NULL
          }

        if (selection.K=="Lav"){
          res.LoopK=loop.iterative.procedure(Data.X,lyear,lmin,Kmax,Used.function,threshold,tol)
          loglik=res.LoopK$loglik
          Kh=MLcriterion(res.LoopK$SSwg, Kseq,S)
          res.segfunct=c()
          request=paste(paste0("res.segfunct=",Used.function,'(Data.X,Kh,lmin,lyear,threshold,tol)'),sep="")
          eval(parse(text=request))
          Tmu=c()
          fh=c()
          varh=c()
          res.segfunct.with.NA<- add_NA(res.segfunct,n.Data,present.data)
          Tmu=res.segfunct.with.NA$Tmu
          fh=res.segfunct.with.NA$f
          varh=res.segfunct$sigma.est.month^2
          mBIC=NULL
        }

        if (selection.K=="BM_BJ"){
          res.LoopK=loop.iterative.procedure(Data.X,lyear,lmin,Kmax,Used.function,threshold,tol)
          loglik=res.LoopK$loglik
          pen=5*Kseq+2*Kseq*log(n.X/Kseq)
          Kh=BMcriterion(res.LoopK$SSwg,pen)
          res.segfunct=c()
          request=paste(paste0("res.segfunct=",Used.function,'(Data.X,Kh,lmin,lyear,threshold,tol)'),sep="")
          eval(parse(text=request))
          Tmu=c()
          fh=c()
          varh=c()
          res.segfunct.with.NA<- add_NA(res.segfunct,n.Data,present.data)
          Tmu=res.segfunct.with.NA$Tmu
          fh=res.segfunct.with.NA$f
          varh=res.segfunct$sigma.est.month^2
          mBIC=NULL
        }

        if (selection.K=="BM_slope"){
          res.LoopK=loop.iterative.procedure(Data.X,lyear,lmin,Kmax,Used.function,threshold,tol)
          loglik=res.LoopK$loglik
          pen=5*Kseq+2*Kseq*log(n.X/Kseq)
          KK=Kmax
          DataForCa=data.frame(model=paste("K=",Kseq[1:KK]),pen=pen[1:KK],complexity=Kseq[1:KK],contrast=res.LoopK$SSwg[1:KK])
          Kh=Kseq[which(capushe::DDSE(DataForCa)@model==DataForCa$model)]
          res.segfunct=c()
          request=paste(paste0("res.segfunct=",Used.function,'(Data.X,Kh,lmin,lyear,threshold,tol)'),sep="")
          eval(parse(text=request))
          Tmu=c()
          fh=c()
          varh=c()
          res.segfunct.with.NA<- add_NA(res.segfunct,n.Data,present.data)
          Tmu=res.segfunct.with.NA$Tmu
          fh=res.segfunct.with.NA$f
          varh=res.segfunct$sigma.est.month^2
          mBIC=NULL
        }

        if (selection.K=="mBIC"){
          res.LoopK=loop.iterative.procedure(Data.X,lyear,lmin,Kmax,Used.function,threshold,tol)
          loglik=res.LoopK$loglik
          bic=mBICcriterion(res.LoopK$SSwg,res.LoopK$LogLg,n.X,Kseq)
          Kh=bic$Kh
          res.segfunct=c()
          request=paste(paste0("res.segfunct=",Used.function,'(Data.X,Kh,lmin,lyear,threshold,tol)'),sep="")
          eval(parse(text=request))
          Tmu=c()
          fh=c()
          varh=c()
          res.segfunct.with.NA<- add_NA(res.segfunct,n.Data,present.data)
          Tmu=res.segfunct.with.NA$Tmu
          fh=res.segfunct.with.NA$f
          varh=res.segfunct$sigma.est.month^2
          mBIC=bic$mBIC
        }

        if (selection.K=="All"){
          res.LoopK=loop.iterative.procedure(Data.X,lyear,lmin,Kmax,Used.function,threshold,tol)
          loglik=res.LoopK$loglik
          Tmu=list()
          Kh=list()
          fh=list()
          varh=list()
          mBIC=c()

           #1=mBIC
          bic=mBICcriterion(res.LoopK$SSwg,res.LoopK$LogLg,n.X,Kseq)
          Kh.mBIC=bic$Kh
          res.segfunct=c()
          request=paste(paste0("res.segfunct=",Used.function,'(Data.X,Kh.mBIC,lmin,lyear,threshold,tol)'),sep="")
          eval(parse(text=request))
          res.segfunct.mBIC<- add_NA(res.segfunct,n.Data,present.data)
          varh=res.segfunct$sigma.est.month^2
          Tmu$mBIC=res.segfunct.mBIC$Tmu
          fh$mBIC=res.segfunct.mBIC$f
          Kh$mBIC=Kh.mBIC
          mBIC=bic$mBIC

          #2=ML
          Kh.ML=MLcriterion(res.LoopK$SSwg, Kseq,S)
          Kh$Lav=Kh.ML
          if (Kh.ML==Kh.mBIC){
            res.segfunct.ML=res.segfunct.mBIC
          } else{
            res.segfunct=c()
            request=paste(paste0("res.segfunct=",Used.function,'(Data.X,Kh.ML,lmin,lyear,threshold,tol)'),sep="")
            eval(parse(text=request))
            res.segfunct.ML<- add_NA(res.segfunct,n.Data,present.data)
          }
          Tmu$Lav=res.segfunct.ML$Tmu
          fh$Lav=res.segfunct.ML$f

          #3=BM_BJ
          pen=5*Kseq+2*Kseq*log(n.X/Kseq)
          Kh.BM_BJ=BMcriterion(res.LoopK$SSwg,pen)
          Kh$BM_BJ=Kh.BM_BJ
          if ((Kh.BM_BJ==Kh.mBIC)) {
            res.segfunct.BM_BJ=res.segfunct.mBIC
          } else if ((Kh.BM_BJ==Kh.ML)) {
            res.segfunct.BM_BJ=res.segfunct.ML
          } else {
            res.segfunct=c()
            request=paste(paste0("res.segfunct=",Used.function,'(Data.X,Kh.BM_BJ,lmin,lyear,threshold,tol)'),sep="")
            eval(parse(text=request))
            res.segfunct.BM_BJ<- add_NA(res.segfunct,n.Data,present.data)
          }
          Tmu$BM_BJ=res.segfunct.BM_BJ$Tmu
          fh$BM_BJ=res.segfunct.BM_BJ$f

          #4=BM2
          KK=Kmax
          DataForCa=data.frame(model=paste("K=",Kseq[1:KK]),pen=pen[1:KK],complexity=Kseq[1:KK],contrast=res.LoopK$SSwg[1:KK])
          Kh.BM_slope=Kseq[which(capushe::DDSE(DataForCa)@model==DataForCa$model)]
          Kh$BM_slope=Kh.BM_slope
          res.segfunct.BM_slope=c()

          if ((Kh.BM_slope==Kh.mBIC)) {
            res.segfunct.BM_slope=res.segfunct.mBIC
          } else if ((Kh.BM_slope==Kh.ML)) {
            res.segfunct.BM_slope=res.segfunct.ML
          } else if ((Kh.BM_slope==Kh.BM_BJ)) {
            res.segfunct.BM_slope=res.segfunct.BM_BJ
          } else {
            res.segfunct=c()
            request=paste(paste0("res.segfunct=",Used.function,'(Data.X,Kh.BM_slope,lmin,lyear,threshold,tol)'),sep="")
            eval(parse(text=request))
            res.segfunct.BM_slope<- add_NA(res.segfunct,n.Data,present.data)
          }
          Tmu$BM_slope=res.segfunct.BM_slope$Tmu
          fh$BM_slope=res.segfunct.BM_slope$f
          }

        result$selected.K=Kh
        result$segmentation=Tmu
        result$functional=fh
        result$loglik=loglik
        result$variances=varh
        result$mBIC=mBIC
        return(result)
        } else {

          sigma.est.month=RobEstiMonthlyVariance(Data.X)
          var.est.k=c()
          var.est.k=sigma.est.month^2
          var.est.t=c()
          var.est.t=var.est.k[as.numeric(Data.X$month)]



          if (selection.K=="none"){
            res.seg=c()
            res.seg=SegMonthlyVarianceK(Data.X,Kmax,lmin,var.est.t)
            res.seg$f=rep(0,n.X)
            varh=sigma.est.month^2
            loglik=-(n.X/2)*(log(2*pi))-(1/2)*(sum(log(var.est.t)))-(1/2)*res.seg$SSwg
            Tmu=c()
            fh=c()
            res.seg.with.NA<- add_NA(res.seg,n.Data,present.data)
            Tmu=res.seg.with.NA$Tmu
            fh=res.seg.with.NA$f
            Kh=Kmax
            mBIC=NULL
          }

          if (selection.K=="Lav"){
            res.seg=c()
            res.seg=SegMonthlyVarianceK(Data.X,Kmax,lmin,var.est.t)
            varh=sigma.est.month^2
            loglik=-(n.X/2)*(log(2*pi))-(1/2)*(sum(log(var.est.t)))-(1/2)*res.seg$SSwg
            Kh=MLcriterion(res.seg$SSwg, Kseq,S)
            res.seg=c()
            res.seg=SegMonthlyVarianceK(Data.X,Kh,lmin,var.est.t)
            res.seg$f=rep(0,n.X)
            Tmu=c()
            fh=c()
            res.seg.with.NA<- add_NA(res.seg,n.Data,present.data)
            Tmu=res.seg.with.NA$Tmu
            fh=res.seg.with.NA$f
            mBIC=NULL
          }

          if (selection.K=="BM_BJ"){
            res.seg=c()
            res.seg=SegMonthlyVarianceK(Data.X,Kmax,lmin,var.est.t)
            varh=sigma.est.month^2
            loglik=-(n.X/2)*(log(2*pi))-(1/2)*(sum(log(var.est.t)))-(1/2)*res.seg$SSwg
            pen=5*Kseq+2*Kseq*log(n.X/Kseq)
            Kh=BMcriterion(res.seg$SSwg,pen)
            res.seg=c()
            res.seg=SegMonthlyVarianceK(Data.X,Kh,lmin,var.est.t)
            res.seg$f=rep(0,n.X)
            Tmu=c()
            fh=c()
            res.seg.with.NA<- add_NA(res.seg,n.Data,present.data)
            Tmu=res.seg.with.NA$Tmu
            fh=res.seg.with.NA$f
            mBIC=NULL
          }

          if (selection.K=="BM_slope"){
            res.seg=c()
            res.seg=SegMonthlyVarianceK(Data.X,Kmax,lmin,var.est.t)
            varh=sigma.est.month^2
            loglik=-(n.X/2)*(log(2*pi))-(1/2)*(sum(log(var.est.t)))-(1/2)*res.seg$SSwg
            pen=5*Kseq+2*Kseq*log(n.X/Kseq)
            KK=Kmax
            DataForCa=data.frame(model=paste("K=",Kseq[1:KK]),pen=pen[1:KK],complexity=Kseq[1:KK],contrast=res.seg$SSwg[1:KK])
            Kh=Kseq[which(capushe::DDSE(DataForCa)@model==DataForCa$model)]
            res.seg=c()
            res.seg=SegMonthlyVarianceK(Data.X,Kh,lmin,var.est.t)
            res.seg$f=rep(0,n.X)
            Tmu=c()
            fh=c()
            res.seg.with.NA<- add_NA(res.seg,n.Data,present.data)
            Tmu=res.seg.with.NA$Tmu
            fh=res.seg.with.NA$f
            mBIC=NULL
          }

          if (selection.K=="mBIC"){
            res.seg=c()
            res.seg=SegMonthlyVarianceK(Data.X,Kmax,lmin,var.est.t)
            varh=sigma.est.month^2
            loglik=-(n.X/2)*(log(2*pi))-(1/2)*(sum(log(var.est.t)))-(1/2)*res.seg$SSwg
            LogLg=c()
            for (k in 1:Kmax){
            seg       = matrix(res.seg$breaks[k,1:k],ncol=k,nrow=1)
            LogLg[k] =  apply(seg,1,FUN=function(z) sum(log(diff(c(0,z))[diff(c(0,z))>0])))
            }
            bic=mBICcriterion(res.seg$SSwg,LogLg,n.X,Kseq)
            Kh=bic$Kh
            res.seg=c()
            res.seg=SegMonthlyVarianceK(Data.X,Kh,lmin,var.est.t)
            res.seg$f=rep(0,n.X)
            Tmu=c()
            fh=c()
            res.seg.with.NA<- add_NA(res.seg,n.Data,present.data)
            Tmu=res.seg.with.NA$Tmu
            fh=res.seg.with.NA$f
            mBIC=bic$mBIC
          }
          if (selection.K=="All"){
            Tmu=list()
            Kh=list()
            fh=list()
            varh=list()
            mBIC=c()
            res.seg=c()
            res.seg=SegMonthlyVarianceK(Data.X,Kmax,lmin,var.est.t)
            varh=sigma.est.month^2
            loglik=-(n.X/2)*(log(2*pi))-(1/2)*(sum(log(var.est.t)))-(1/2)*res.seg$SSwg


            #1=mBIC
            LogLg=c()
            for (k in 1:Kmax){
              seg       = matrix(res.seg$breaks[k,1:k],ncol=k,nrow=1)
              LogLg[k] =  apply(seg,1,FUN=function(z) sum(log(diff(c(0,z))[diff(c(0,z))>0])))
            }
            bic=mBICcriterion(res.seg$SSwg,LogLg,n.X,Kseq)
            Kh.mBIC=bic$Kh
            res.seg.temp=c()
            res.seg.temp=SegMonthlyVarianceK(Data.X,Kh.mBIC,lmin,var.est.t)
            res.seg.temp$f=rep(0,n.X)
            res.seg.mBIC<- add_NA(res.seg.temp,n.Data,present.data)
            Tmu$mBIC=res.seg.mBIC$Tmu
            fh$mBIC=res.seg.mBIC$f
            Kh$mBIC=Kh.mBIC
            mBIC=bic$mBIC

            #2=ML
            Kh.ML=MLcriterion(res.seg$SSwg, Kseq,S)
            Kh$Lav=Kh.ML
            if (Kh.ML==Kh.mBIC){
              res.seg.ML=res.seg.mBIC
            } else{
              res.seg.temp=c()
              res.seg.temp=SegMonthlyVarianceK(Data.X,Kh.ML,lmin,var.est.t)
              res.seg.temp$f=rep(0,n.X)
              res.seg.ML<- add_NA(res.seg.temp,n.Data,present.data)
            }
            Tmu$Lav=res.seg.ML$Tmu
            fh$Lav=res.seg.ML$f

            #3=BM_BJ
            pen=5*Kseq+2*Kseq*log(n.X/Kseq)
            Kh.BM_BJ=BMcriterion(res.seg$SSwg,pen)
            Kh$BM_BJ=Kh.BM_BJ
            if ((Kh.BM_BJ==Kh.mBIC)) {
              res.seg.BM_BJ=res.seg.mBIC
            } else if ((Kh.BM_BJ==Kh.ML)) {
              res.seg.BM_BJ=res.seg.ML
            } else {
              res.seg.temp=c()
              res.seg.temp=SegMonthlyVarianceK(Data.X,Kh.BM_BJ,lmin,var.est.t)
              res.seg.temp$f=rep(0,n.X)
              res.seg.BM_BJ<- add_NA(res.seg.temp,n.Data,present.data)
            }
            Tmu$BM_BJ=res.seg.BM_BJ$Tmu
            fh$BM_BJ=res.seg.BM_BJ$f

            #4=BM2
            KK=Kmax
            DataForCa=data.frame(model=paste("K=",Kseq[1:KK]),pen=pen[1:KK],complexity=Kseq[1:KK],contrast=res.seg$SSwg[1:KK])
            Kh.BM_slope=Kseq[which(capushe::DDSE(DataForCa)@model==DataForCa$model)]
            Kh$BM_slope=Kh.BM_slope
            res.seg.BM_slope=c()

            if ((Kh.BM_slope==Kh.mBIC)) {
              res.seg.BM_slope=res.seg.mBIC
            } else if ((Kh.BM_slope==Kh.ML)) {
              res.seg.BM_slope=res.seg.ML
            } else if ((Kh.BM_slope==Kh.BM_BJ)) {
              res.seg.BM_slope=res.seg.BM_BJ
            } else {
              res.seg.temp=c()
              res.seg.temp=SegMonthlyVarianceK(Data.X,Kh.BM_slope,lmin,var.est.t)
              res.seg.temp$f=rep(0,n.X)
              res.segKh.BM_slope<- add_NA(res.seg.temp,n.Data,present.data)
            }
            Tmu$BM_slope=res.seg.BM_slope$Tmu
            fh$BM_slope=res.seg.BM_slope$f
          }

          result$selected.K=Kh
          result$segmentation=Tmu
          result$functional=fh
          result$loglik=loglik
          result$variances=varh
          result$mBIC=mBIC
          return(result)
        }
  }
  }


