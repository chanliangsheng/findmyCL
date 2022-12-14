% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rtFilter.R
\name{delete_no_FA_MS1}
\alias{delete_no_FA_MS1}
\title{Delete redundant Cardiolipin.}
\usage{
delete_no_FA_MS1(bigMS1)
}
\arguments{
\item{bigMS1}{list(1),a bigMS1 contain MS1,MS2,splice-Cardiolipin-result.}
}
\value{
list(1),a bigMS1 contain MS1,MS2,splice-Cardiolipin-result.
}
\description{
Delete redundant Cardiolipin and splicePA-result in MS1 which can not spliced by PA&FA.
}
