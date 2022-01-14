loadTestData <- function(){
  rawdata <- findmyCL::loadCentroidData(file = "50X_NEG_CID_150-2000_3uL(centroid).mzml") %>%
    findmyCL::matchMS1(ppm = 10) %>%
    findmyCL::checkMS2(ppm = 5) %>%
    return()
}
