% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matchMS1.R
\name{matchMS1}
\alias{matchMS1}
\title{Matching the MS1 database}
\usage{
matchMS1(object, ppm = 5, database = c("CL", "MLCL", "DLCL"))
}
\arguments{
\item{object}{a findCL object from \code{\link[=loadCentroidData]{loadCentroidData()}}}

\item{ppm}{numeric(1) defining the maximal tolerated m/z deviation in consecutive scans in parts per million (ppm) for the initial ROI definition when match our database.Default is 5.}

\item{database}{vector,choose a database you like to match your data,like("CL","MLCL","DLCL")}
}
\value{
a findmyCL object
}
\description{
use database to select which MS1 is more likely to be a Cardiolipin.
}
\examples{
setwd("D:/mzml file")
# set the path of your mzml file

rawdata <- findmyCL::loadCentroidData(file = "50X_NEG_CID_150-2000_3uL(centroid).mzml" , ppm = 5)
# read the mzml file in centroid mode

MS1_match_result <- findmyCL::matchMS1(object = rawdata , ppm = 5, database = c("CL","MLCL","DLCL"))
# matching MS1 database,we use all database here.
}
\seealso{
\code{\link[=loadCentroidData]{loadCentroidData()}} create a findCL object to be an input for this function
}
