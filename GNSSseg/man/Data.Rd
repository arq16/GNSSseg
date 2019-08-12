\name{Data}
\alias{Data}
\docType{data}
\title{Example of data}
\usage{
data(Data)
}
\description{A data frame [n x 2] containing a simulated Gaussian series for the two years 1995 and 1996, with size n=731. 3 changes are considered at positions 100, 150 and 500 (or at the dates 1995-04-10, 1995-05-30 and 1996-05-14). The means of the segments alternates between 0 and 1 (beginning by 0). The functional part is 0.4*cos(2*pi*time/lyear) where lyear is 365.25 and time is centred according to the first date and expressed in seconds: time=(date-date[1])/86400. The standard deviation of the noise of the 12 months are drawn from an uniform distribution between 0.1 and 0.8. The date is expressed as yyyy-mm-dd in the "calendar time" (class POSIXct).
}
\format{
  A data frame with 731 observations on the following 2 variables.
  \describe{
    \item{\code{signal}}{a numeric vector}
    \item{\code{date}}{a date vector expressed as yyyy-mm-dd in the "calendar time" (class POSIXct)}
  }
}
\examples{
library(GNSSseg)
data(Data)
class(Data$date)
plot(Data$date,Data$signal,type="l")
}
\details{signal: the values of the observed signal; date: the dates in calendar time
}
\keyword{datasets}
