#' @title Data preprocessing
#' @description  A function to protreat the mzml file in profile mode into a mzml file in centroid mode,the result of this function is to output a mzml file in centroid mode.https://bioconductor.riken.jp/packages/3.7/bioc/vignettes/MSnbase/inst/doc/MSnbase-centroiding.html
#' @importFrom MSnbase readMSData smooth pickPeaks writeMSData
#' @param file character(1)  mzml file in profile mode
#' @param smooth_method character(1) smoothing method for your data.Like MovingAverage and SavitzkyGolay.
#' @param halfWindowSize numeric(1) if you choose MovingAverage for your smooth,you have to set the halfWindowSize,default is 2.
#' @param polynomialOrder character(1) if you choose SavitzkyGolay for your smooth,you have to set the polynomialOrder,default is 3.
#' @param pickPeaks_method character(1) peak picking method for your data.Like kNeighbors and descendPeak.
#' @param kNeighbors_par numeric(1) if you choose kNeighbors for your peak picking method,you have to set the parameter,kNeighbors_par,default is 1.
#' @param signalPercentage numeric(1) if you choose descendPeak for your peak picking method,you have to set the parameter,signalPercentage,default is 50.
#' @param output_file character(1) the name of your output mzml file.
#' @return a mzml file
#' @examples
#' profilePretreat(file = "QC-NGE-ALL.mzML" , smooth_method = "MovingAverage" , halfWindowSize = 2 , pickPeaks_method = "kNeighbors" , kNeighbors_par = 1 , output_file = "Test.mzml")
#' @export
profilePretreat <- function(
  file,
  smooth_method = "MovingAverage",
  halfWindowSize = 2,
  polynomialOrder = 3,
  pickPeaks_method = "kNeighbors",
  kNeighbors_par = 1,
  signalPercentage = 50,
  output_file){

  message("Reaing...")
  data_prof_inmemory <- readMSData(file, mode = "onDisk", centroided = FALSE , smoothed. = FALSE)
  message("Done!")
  #读取mzml文件，模式为inMemory，先只提取一级峰

  if (smooth_method == "MovingAverage") {
    message("Smoothing...")
    data_prof_smooth <- smooth(data_prof_inmemory, method = smooth_method, halfWindowSize = halfWindowSize)
    #对质谱数据进行平滑，过滤噪声
    message("Done!")
  }
  if (smooth_method == "SavitzkyGolay") {
    message("Smoothing...")
    data_prof_smooth <- smooth(data_prof_inmemory, method = smooth_method, polynomialOrder = polynomialOrder)
    #对质谱数据进行平滑，过滤噪声
    message("Done!")
  }

  if (pickPeaks_method == "kNeighbors") {
    message("pickPeaking...")
    data_prof_pick <- pickPeaks(data_prof_smooth , refinMz = pickPeaks_method, k = kNeighbors_par)
    #对质谱数据进行峰拾取，使profile变为centroid
    #pickpeaks会减少原本数据的大量内存，因为减少了很大的数据量
    message("Done!")
  }
  if (pickPeaks_method == "descendPeak") {
    message("pickPeaking...")
    data_prof_pick <- pickPeaks(data_prof_smooth , refinMz = pickPeaks_method, signalPercentage = signalPercentage)
    #对质谱数据进行峰拾取，使profile变为centroid
    message("Done!")
  }

  paste0("Writing data to ",output_file," ...") %>% message()
  writeMSData(object = data_prof_pick , file = output_file , outformat = "mzml", merge = FALSE , copy = FALSE , software_processing = NULL)
  #将centroid的文件输出为mzml文件
  message("Done!")

}
