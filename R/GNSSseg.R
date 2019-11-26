#' Homogeneization of GNSS series
#'
#' fit a segmentation in the mean model by taken into account for a functional part and a heterogeneous variance (default is monthly)
#'
#' @param Data a data frame, with size [n x 2], containing the signal (e.g. the daily GPS-ERAI series for GNSS) and the dates (in format yyyy-mm-dd of type "calendar time" (class POSIXct))
#' @param lyear the length of the year in the signal. Default is 365.25
#' @param lmin the minimum length of the segments. Default is 1
#' @param Kmax the maximal number of segments (must be lower than n). Default is 30
#' @param selection.K a name indicating the model selection criterion to select the number of segments K (\code{mBIC}, \code{Lav}, \code{BM_BJ} or \code{BM_slope}). \code{"none"} indicates that no selection is claimed and the procedure considers \code{Kmax} segments or \code{Kmax}-1 changes. If \code{selection.K="All"}, the results for the four possible criteria are given. Default is \code{"mBIC"}
#' @param S the threshold used in the Lav's criterion. Default is 0.75
#' @param f a boolean indicating if the functional part is taking into account in the model. Default is TRUE and note that if \code{f=FALSE}, only a segmentation is performed
#' @param selection.f a boolean indicating if a selection on the functions of the Fourier decomposition of order 4 is performed. Default is FALSE
#' @param threshold a numeric value lower than 1 used for the selection of the functions of the Fourier decomposition of order 4. Default is 0.001
#' @param tol the stopping rule for the iterative procedure. Default is 1e-4
#'
#' @return A file containing
#' \itemize{
#' \item \code{K} that corresponds to the selected number of segments or \code{K}-1 corresponds to the number of changes. If \code{selection.K="none"}, the number of segments is \code{Kmax}. 
#' \item \code{seg} that corresponds to the estimation of the segmentation parameters (the begin and the end positions of each segment with the estimated mean). 
#' \item \code{funct} that corresponds to the estimation of the functional part. If \code{f==FALSE}, \code{funct} is FALSE
#' \item \code{coeff} that corresponds to the estimation of the coefficients of the Fourier decomposition. The vector contains 8 coefficients if \code{selection.f=FALSE} or as many coefficients as the number of selected functions if \code{selection.f=TRUE}. If \code{f==FALSE}, \code{coeff} is FALSE
#' \item \code{variances} that corresponds to the estimated variances of each fixed interval
#' \item \code{SSR} that corresponds to the Residuals Sum of Squares for k=1,...,\code{Kmax}. If \code{selection.K="none"}, it contains only the SSR for \code{Kmax} segments
#' \item \code{Tot} is a list. Each component contains all the results k segments (k=1,...,\code{Kmax}). If \code{selection.K="none"}, \code{Tot} is NA
#' }
#' If \code{selection.K="All"}, the outputs \code{K}, \code{seg}, \code{funct} and \code{coeff} are each a list containing the corresponding results obtained for the four model selection criteria
#'
#' @details
#' The function performs homogeneization of GNSS series. The considered model is such that: (1) the average is composed of a piecewise function (changes in the mean) with a functional part and (2) the variance is heterogeneous on fixed intervals. By default the latter intervals are the months. 
#' The inference procedure consists in two steps. First, the number of segments is fixed to \code{Kmax} and the parameters are estimated using the maximum likelihood procedure using the following procedure: first the variances are robustly estimated and then the segmentation and the functional parts are iteratively estimated. Then the number of segments is chosen using model selection criteria. The possible criteria are \code{mBIC} the modified BIC criterion REFEREF, \code{Lav} the criterion proposed by REFEF, \code{BM_BJ} and \code{BM_slope} the criteria proposed by REFEF where the penalty constant is calibrated using the Biggest Jump and the slope respectively REFERF.
#' \itemize{
#' \item The data is a data frame with 2 columns: $signal is the signal to be homogeneized (a daily series) and $date is the date. The date will be in format yyyy-mm-dd of type "calendar time" (class POSIXct).
#' \item The function part is estimated using a Fourier decomposition of order 4 with \code{selection.f=FALSE}. \code{selection.f=TRUE} consists in selecting the significative functions of the Fourier decomposition of order 4 (for which p.values are lower than \code{threshold})
#' \item If \code{selection.K="none"}, the procedure is performed with \code{Kmax} segments.
#' \item Missing data in the signal are accepted.
#' }
#' 
#' @examples
#' data(Data)
#' lyear=365.25
#' Kmax=10
#' lmin=1
#' result=GNSSseg(Data,lyear,Kmax=Kmax,selection.K="none")
#' plot_GNSS(Data,result$seg,result$funct)
#' @export

GNSSseg=function(Data,lyear=365.25,lmin=1,Kmax=30,selection.K="BM_BJ",S=0.75,f=TRUE,selection.f=FALSE,threshold=0.001,tol=1e-4){
  result  = list()
  Data.X  = c()
  cond1=TRUE
  cond2=TRUE
  
  #For NA
  present.data = which(!is.na(Data$signal))
  Data.X       = Data[present.data,]
  n.Data=length(Data$signal)
  n.X=length(Data.X$signal)
  Kseq=1:Kmax
  
  
  #The conditions to be fulfilled
  if (class(Data$date)[1]!="POSIXct"){
    cond1=FALSE
    cat("date must be in format yyyy-mm-dd of type GMT in class POSIXct/POSIXt")
    }
  if (Kmax >n.X) {
    cond2=FALSE
    cat("The maximal number of segments Kmax", Kmax," needs to be lower than the length of the series without NA that is " ,n.present,"\n")
    }
  
  
  if ((cond1==TRUE) & (cond2==TRUE)){
      Data.X$year=as.factor(format(Data.X$date,format='%Y'))
      Data.X$month=as.factor(format(Data.X$date,format='%m'))
      Data.X$month = droplevels(Data.X$month)
      Data.X$year  = droplevels(Data.X$year)
      
      
      #Used function for NA
      add_NA=function(res,present.data,n.Data,segf){
        res.with.NA=list()
        #Segmentation
        Tmu.temp=res$Tmu
        Tmu.temp$begin=present.data[Tmu.temp$begin]
        Tmu.temp$end=present.data[Tmu.temp$end]
        Tmu.temp$end[length(Tmu.temp$end)]=n.Data
        #Function
        if (segf==TRUE){
          f.temp=rep(NA,n.Data)
          f.temp[present.data]=res$f
          res.with.NA$f=f.temp
        } else {f.temp=FALSE}
        res.with.NA$Tmu=Tmu.temp
        return(res.with.NA)
      }


      #Estimation of the Montly variances
      sigma.est.month=RobEstiMonthlyVariance(Data.X)
      var.est.month=sigma.est.month^2
    
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
          request=paste(paste0("res.segfunct=",Used.function,'(Data.X,var.est.month,Kmax,lmin,lyear,threshold,tol)'),sep="")
          eval(parse(text=request))
          res.segfunct.with.NA= add_NA(res.segfunct,present.data,n.Data,segf=TRUE)
          Tmu=res.segfunct.with.NA$Tmu
          funct=res.segfunct.with.NA$f
          Kh=Kmax
          res.LoopK=NA
          coeff=res.segfunct$coeff
          SSwg=res.segfunct$SSwg
          } 

        if (selection.K=="Lav"){
          res.LoopK=Loop.K.procedure(Data.X,var.est.month,lyear,lmin,Kmax,Used.function,threshold,tol)
          res=sapply(res.LoopK,function(e) {
            return(c(SSwg =e$SSwg, LogLg = e$LogLg))
          })
          SSwg=res[1,]
          Kh=MLcriterion(SSwg, Kseq,S)
          res.segfunct=res.LoopK[[Kh]]
          res.segfunct.with.NA<- add_NA(res.segfunct,present.data,n.Data,segf=TRUE)
          Tmu=res.segfunct.with.NA$Tmu
          funct=res.segfunct.with.NA$f
          coeff=res.segfunct$coeff
        }

        if (selection.K=="BM_BJ"){
          res.LoopK=Loop.K.procedure(Data.X,var.est.month,lyear,lmin,Kmax,Used.function,threshold,tol)
          res=sapply(res.LoopK,function(e) {
            return(c(SSwg =e$SSwg, LogLg = e$LogLg))
          })
          SSwg=res[1,]
          pen=5*Kseq+2*Kseq*log(n.X/Kseq)
          Kh=BMcriterion(SSwg,pen)
          res.segfunct=res.LoopK[[Kh]]
          res.segfunct.with.NA<- add_NA(res.segfunct,present.data,n.Data,segf=TRUE)
          Tmu=res.segfunct.with.NA$Tmu
          funct=res.segfunct.with.NA$f
          coeff=res.segfunct$coeff
        }

        if (selection.K=="BM_slope"){
          res.LoopK=Loop.K.procedure(Data.X,var.est.month,lyear,lmin,Kmax,Used.function,threshold,tol)
          res=sapply(res.LoopK,function(e) {
            return(c(SSwg =e$SSwg, LogLg = e$LogLg))
          })
          SSwg=res[1,]
          pen=5*Kseq+2*Kseq*log(n.X/Kseq)
          DataForCa=data.frame(model=paste("K=",Kseq),pen=pen,complexity=Kseq,contrast=SSwg)
          Kh=Kseq[which(capushe::DDSE(DataForCa)@model==DataForCa$model)]
          res.segfunct=res.LoopK[[Kh]]
          res.segfunct.with.NA<- add_NA(res.segfunct,present.data,n.Data,segf=TRUE)
          Tmu=res.segfunct.with.NA$Tmu
          funct=res.segfunct.with.NA$f
          coeff=res.segfunct$coeff
        }

        if (selection.K=="mBIC"){
          res.LoopK=Loop.K.procedure(Data.X,var.est.month,lyear,lmin,Kmax,Used.function,threshold,tol)
          res=sapply(res.LoopK,function(e) {
            return(c(SSwg =e$SSwg, LogLg = e$LogLg))
          })
          SSwg=res[1,]
          LogLg=res[2,]
          Kh=mBICcriterion(SSwg,LogLg,n.X,Kseq)$Kh
          res.segfunct=res.LoopK[[Kh]]
          res.segfunct.with.NA<- add_NA(res.segfunct,present.data,n.Data,segf=TRUE)
          Tmu=res.segfunct.with.NA$Tmu
          funct=res.segfunct.with.NA$f
          coeff=res.segfunct$coeff
        }

        if (selection.K=="All"){
          res.LoopK=Loop.K.procedure(Data.X,var.est.month,lyear,lmin,Kmax,Used.function,threshold,tol)
          res=sapply(res.LoopK,function(e) {
            return(c(SSwg =e$SSwg, LogLg = e$LogLg))
          })
          SSwg=res[1,]
          LogLg=res[2,]
          pen=5*Kseq+2*Kseq*log(n.X/Kseq)
          Tmu=list()
          Kh=list()
          funct=list()
          coeff=list()

           #1=mBIC
          Kh$mBIC=mBICcriterion(SSwg,LogLg,n.X,Kseq)$Kh
          res.segfunct=c()
          res.segfunct=res.LoopK[[Kh$mBIC]]
          res.segfunct.mBIC= add_NA(res.segfunct,present.data,n.Data,segf=TRUE)
          Tmu$mBIC=res.segfunct.mBIC$Tmu
          funct$mBIC=res.segfunct.mBIC$f
          coeff$mBIC=res.segfunct$coeff
          
          #2=ML
          Kh$Lav=MLcriterion(SSwg, Kseq,S)
          res.segfunct=c()
          res.segfunct=res.LoopK[[Kh$Lav]]
          res.segfunct.ML= add_NA(res.segfunct,present.data,n.Data,segf=TRUE)
          Tmu$Lav=res.segfunct.ML$Tmu
          funct$Lav=res.segfunct.ML$f
          coeff$Lav=res.segfunct$coeff

          #3=BM_BJ
          Kh$BM_BJ=BMcriterion(SSwg,pen)
          res.segfunct=c()
          res.segfunct=res.LoopK[[Kh$BM_BJ]]
          res.segfunct.BM_BJ= add_NA(res.segfunct,present.data,n.Data,segf=TRUE)
          Tmu$BM_BJ=res.segfunct.BM_BJ$Tmu
          funct$BM_BJ=res.segfunct.BM_BJ$f
          coeff$BM_BJ=res.segfunct$coeff
          
          #4=BM2
          DataForCa=data.frame(model=paste("K=",Kseq),pen=pen,complexity=Kseq,contrast=SSwg)
          Kh$BM_slope=Kseq[which(capushe::DDSE(DataForCa)@model==DataForCa$model)]
          res.segfunct=c()
          res.segfunct=res.LoopK[[Kh$BM_slope]]
          res.segfunct.BM_slope= add_NA(res.segfunct,present.data,n.Data,segf=TRUE)
          Tmu$BM_slope=res.segfunct.BM_slope$Tmu
          funct$BM_slope=res.segfunct.BM_slope$f
          coeff$BM_slope=res.segfunct$coeff
          
          }

        
        
        } else {
          funct=FALSE
          coeff=FALSE
          var.est.t=var.est.month[as.numeric(Data.X$month)]
          res.seg=SegMonthlyVarianceK(Data.X,Kmax,lmin,var.est.t)
          SSwg=res.seg$SSwg
          pen=5*Kseq+2*Kseq*log(n.X/Kseq)
          res.LoopK=res.seg$res.LoopK
          
          if (selection.K=="none"){
            res.seg.with.NA<- add_NA(res.seg,present.data,n.Data,segf=FALSE)
            Tmu=res.seg.with.NA$Tmu
            Kh=Kmax
          }

          if (selection.K=="Lav"){
            Kh=MLcriterion(SSwg, Kseq,S)
            res.seg.sol=c()
            res.seg.sol=SegMonthlyVarianceK(Data.X,Kh,lmin,var.est.t)
            res.seg.with.NA<- add_NA(res.seg.sol,present.data,n.Data,segf=FALSE)
            Tmu=res.seg.with.NA$Tmu
          }

          if (selection.K=="BM_BJ"){
            Kh=BMcriterion(SSwg,pen)
            res.seg.sol=c()
            res.seg.sol=SegMonthlyVarianceK(Data.X,Kh,lmin,var.est.t)
            res.seg.with.NA<- add_NA(res.seg.sol,present.data,n.Data,segf=FALSE)
            Tmu=res.seg.with.NA$Tmu
          }

          if (selection.K=="BM_slope"){
            DataForCa=data.frame(model=paste("K=",Kseq),pen=pen,complexity=Kseq,contrast=SSwg)
            Kh=Kseq[which(capushe::DDSE(DataForCa)@model==DataForCa$model)]
            res.seg.sol=c()
            res.seg.sol=SegMonthlyVarianceK(Data.X,Kh,lmin,var.est.t)
            res.seg.with.NA<- add_NA(res.seg.sol,present.data,n.Data,segf=FALSE)
            Tmu=res.seg.with.NA$Tmu
          }

          if (selection.K=="mBIC"){
            Kh=mBICcriterion(SSwg,res.seg$LogLg,n.X,Kseq)$Kh
            res.seg.sol=c()
            res.seg.sol=SegMonthlyVarianceK(Data.X,Kh,lmin,var.est.t)
            res.seg.with.NA<- add_NA(res.seg.sol,present.data,n.Data,segf=FALSE)
            Tmu=res.seg.with.NA$Tmu
          }
          
          if (selection.K=="All"){
            Tmu=list()
            Kh=list()
            
            #1=mBIC
            Kh.mBIC=mBICcriterion(SSwg,res.seg$LogLg,n.X,Kseq)$Kh
            res.seg.sol=c()
            res.seg.sol=SegMonthlyVarianceK(Data.X,Kh.mBIC,lmin,var.est.t)
            res.seg.mBIC<- add_NA(res.seg.sol,present.data,n.Data,segf=FALSE)
            Tmu$mBIC=res.seg.mBIC$Tmu
            Kh$mBIC=Kh.mBIC
            
            #2=ML
            Kh.ML=MLcriterion(SSwg, Kseq,S)
            Kh$Lav=Kh.ML
            if (Kh.ML==Kh.mBIC){
              res.seg.ML=res.seg.mBIC
            } else{
              res.seg.sol=c()
              res.seg.sol=SegMonthlyVarianceK(Data.X,Kh.ML,lmin,var.est.t)
              res.seg.ML<- add_NA(res.seg.sol,present.data,n.Data,segf=FALSE)
            }
            Tmu$Lav=res.seg.ML$Tmu

            #3=BM_BJ
            Kh.BM_BJ=BMcriterion(SSwg,pen)
            Kh$BM_BJ=Kh.BM_BJ
            if ((Kh.BM_BJ==Kh.mBIC)) {
              res.seg.BM_BJ=res.seg.mBIC
            } else if ((Kh.BM_BJ==Kh.ML)) {
              res.seg.BM_BJ=res.seg.ML
            } else {
              res.seg.sol=c()
              res.seg.sol=SegMonthlyVarianceK(Data.X,Kh.BM_BJ,lmin,var.est.t)
              res.seg.BM_BJ<- add_NA(res.seg.sol,present.data,n.Data,segf=FALSE)
            }
            Tmu$BM_BJ=res.seg.BM_BJ$Tmu

            #4=BM2
            
            DataForCa=data.frame(model=paste("K=",Kseq),pen=pen,complexity=Kseq,contrast=SSwg)
            Kh.BM_slope=Kseq[which(capushe::DDSE(DataForCa)@model==DataForCa$model)]
            Kh$BM_slope=Kh.BM_slope
            if ((Kh.BM_slope==Kh.mBIC)) {
              res.seg.BM_slope=res.seg.mBIC
            } else if ((Kh.BM_slope==Kh.ML)) {
              res.seg.BM_slope=res.seg.ML
            } else if ((Kh.BM_slope==Kh.BM_BJ)) {
              res.seg.BM_slope=res.seg.BM_BJ
            } else {
              res.seg.sol=c()
              res.seg.sol=SegMonthlyVarianceK(Data.X,Kh.BM_slope,lmin,var.est.t)
              res.seg.BM_slope= add_NA(res.seg.sol,present.data,n.Data,segf=FALSE)
            }
            Tmu$BM_slope=res.seg.BM_slope$Tmu
          }
          }
      
      #Obtained segmentation
      result$K=Kh
      result$seg=Tmu
      result$funct=funct
      result$coeff=coeff
      #Global results
      result$variances=var.est.month
      result$SSR=SSwg
      result$Tot=res.LoopK
      return(result)
      }

}

