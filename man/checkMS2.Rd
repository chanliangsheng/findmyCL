% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/checkMS2.R
\name{checkMS2}
\alias{checkMS2}
\title{Checking which MS1 have MS2,from\code{\link[=matchMS1]{matchMS1()}}}
\usage{
checkMS2(object, ppm = 5)
}
\arguments{
\item{object}{a findmyCL object}

\item{ppm}{numeric(1) defining the maximal tolerated m/z deviation in consecutive scans in parts per million (ppm) for the initial ROI definition when matching with the precursorMz of MS2,default is 5.}
}
\value{
a findmyCL object
}
\description{
Checking if a MS1 have corresponding MS2,use this function after \code{\link[=matchMS1]{matchMS1()}}.A MS1 may have more than 1 MS2.The result is stored in ms2CheckResult.
}
