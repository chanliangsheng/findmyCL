#' @title Load rawdata in centroid mode
#' @description A function use function xcmsSet of  package "xcms" and "MSnbase" for deceting the peak of mzml file in centroid mode.We set some defaults in this function like ppm and snthresh.If your data format is not centroid,you can use [findmyCL::profilePretreat()]
#' @importFrom xcms xcmsSet
#' @importFrom MSnbase readMSData
#' @param file a mzml file in centroid mode
#' @param ppm numeric(1) defining the maximal tolerated m/z deviation in consecutive scans in parts per million (ppm) for the initial ROI definition.Default is 5.
#' @param peakwidth numeric(2) with the expected approximate peak width in chromatographic space. Given as a range (min, max) in seconds.Default is c(5,30).
#' @param prefilter numeric(2): c(k, I) specifying the prefilter step for the first analysis step (ROI detection). Mass traces are only retained if they contain at least k peaks with intensity >= I.Default is c(4,5000).
#' @param snthresh  numeric(1) defining the signal to noise ratio cutoff.Default is 3.
#' @param mzdiff numeric(1) representing the minimum difference in m/z dimension required for peaks with overlapping retention times; can be negative to allow overlap. During peak post-processing, peaks defined to be overlapping are reduced to the one peak with the largest signal.Default is -0.001.
#' @param noise numeric(1) allowing to set a minimum intensity required for centroids to be considered in the first analysis step (centroids with intensity < noise are omitted from ROI detection).Default is 0.
#' @param integrate Integration method. For integrate = 1 peak limits are found through descent on the mexican hat filtered data, for integrate = 2 the descent is done on the real data. The latter method is more accurate but prone to noise, while the former is more robust, but less exact.Default is 1.
#' @param fitgauss logical(1) whether or not a Gaussian should be fitted to each peak. This affects mostly the retention time position of the peak.Default is FALSE.
#' @param MS2cutRadio numeric(1) cut the last intensity in MS2
#' @param rtcut numeric(1) half peak width to cut the MS2
#' @seealso [profilePretreat()] turn a mzml file in profile mode to  centroid mode for the input of this function
#' @return a findmyCL object
#' @export
loadCentroidData <- function(file , ppm = 5 , peakwidth = c(5, 30) , prefilter = c(4, 5000) , snthresh = 3 , mzdiff = -0.001 , noise= 0 , integrate = 1L , fitgauss = FALSE , MS2cutRadio = 0.1 , rtcut = 6){

  defineMyClass()
  #定义findmyCL类

  message("Picking peak...")
  rawdata <- xcmsSet(file, method = 'centWave', ppm = ppm,peakwidth = peakwidth,prefilter = prefilter,snthresh = snthresh,mzdiff = mzdiff,noise= noise,mzCenterFun = "wMean",integrate = integrate,fitgauss = fitgauss,mslevel = 1)
  #进行峰拾取
  peak <- rawdata@peaks %>%
    as.data.frame()
  peak$rtmax <- peak$rt + rtcut
  peak$rtmin <- peak$rt - rtcut
  #rt加减定义的半峰宽,这里rtmin可能为负值，不影响后续进展，这里作用是判断哪些二级的保留时间落入一级的保留时间中，如果落入为[负数，正数]中，则必落入[0，正数]中.
  rawdata@peaks <- as.matrix(peak)
  #赋值回原来对象
  message("Done!")
  message("Picking MS2...,please wait...")
  Ms2_data <- MSnbase::readMSData(files = file , msLevel. = 2 ,centroided. = FALSE , smoothed. = FALSE , mode = "inMemory")
  #读取对应文件的二级质谱

  Ms2_data <- findmyCL::cutMS2(MS2 = Ms2_data , MS2cutRadio = MS2cutRadio)
  #将提取的MS2位于最后百分之10的部分去除(可自行设定比率)
  message("Done!")

  findCLobject <- new(Class = "findmyCLclass",path = getwd() , file = file , xcms = rawdata)
  #实例化，生成findmyCL类，并将数据放入里面

  findCLobject@MS2 <- Ms2_data
  #将读取出来的二级质谱数据写入对象中
  findCLobject@loadCentroidData_parameter <- list("ppm" = ppm , "peakwidth" = peakwidth,"snthresh" = snthresh,"mzdiff" = mzdiff , "noise" = noise,"mzCenterFun" = "wMean","integrate" = integrate,"fitgauss" = fitgauss)
  #加入参数的显示

  return(findCLobject)
  #返回这个对象
}

#' @title Delete the peak with the lowest intensity in all MS2.
#' @param MS2 environment(1)
#' @param MS2cutRadio numeric(1),default is 10%.
#' @return environment(1)
#' @seealso [findmyCL::cutMS2_main()]
#' @export
cutMS2 <- function(MS2 , MS2cutRadio = 0.1){
  list <- as.list(MS2@assayData)
  #将变量从环境类型转换为列表类型

  assayData <- purrr::map(.x = list , .f = findmyCL::cutMS2_main , MS2cutRadio = MS2cutRadio) %>%
    as.environment()
  #将最低的百分之10去除

  MS2@assayData <- assayData
  #赋值到MS2中

  return(MS2)
  #返回提取的MS2
}
#将提取的MS2位于最后百分之10的部分去除

#' @title Delete the peak with the lowest intensity in a MS2.
#' @param siMS2 list(1)
#' @param MS2cutRadio numeric(1),default is 10%.
#' @return list(1)
#' @export
cutMS2_main <- function(siMS2 , MS2cutRadio = 0.1){
  delete_length <- length(siMS2@intensity) * MS2cutRadio
  delete_length <- floor(delete_length)
  #求需要删除多少个mz
  if (delete_length == 0) {
    return(siMS2)
  }
  #如果mz总个数小于10，则不删除mz

  delete_number <- 1:delete_length
  #获得排序之后要删除的序列号

  delete <- match(delete_number , order(siMS2@intensity))
  #看intensity里面的哪些是去除的
  delete_result <- order(siMS2@intensity)[-delete]
  #去除最小之后的结果

  siMS2@intensity <- siMS2@intensity[delete_result]
  siMS2@mz <- siMS2@mz[delete_result]
  #结果赋值

  return(siMS2)
}
#将单个二级的位于最后百分之10的部分去除
