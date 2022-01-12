#' @title Matching the MS1 database
#' @description use database to select which MS1 is more likely to be a Cardiolipin.
#' @import furrr
#' @import purrr
#' @import dplyr
#' @importFrom plyr rbind.fill
#' @param object a findCL object from [loadCentroidData()]
#' @param ppm numeric(1) defining the maximal tolerated m/z deviation in consecutive scans in parts per million (ppm) for the initial ROI definition when match our database.Default is 5.
#' @param database vector,choose a database you like to match your data,like("CL","MLCL","DLCL")
#' @seealso [loadCentroidData()] create a findCL object to be an input for this function
#' @return a findCL object
#' @export
#' @examples
#' setwd("D:/mzml file")
#' # set the path of your mzml file
#'
#' rawdata <- findmyCL::loadCentroidData(file = "50X_NEG_CID_150-2000_3uL(centroid).mzml" , ppm = 5)
#' # read the mzml file in centroid mode
#'
#' MS1_match_result <- findmyCL::matchMS1(object = rawdata , ppm = 5, database = c("CL","MLCL","DLCL"))
#' # matching MS1 database,we use all database here.
matchMS1 <- function(object , ppm = 5 , database = c("CL","MLCL","DLCL")){
  database <- as.list(database)
  #转换为列表用于批量操作
  ms1MatchResult <- purrr::map(.x = database,.f = findmyCL::matchMS1_main,object = object , ppm = ppm)
  #计算所有数据配对的结果
  names(ms1MatchResult) <- database
  #重命名结果
  ms1MatchResult <- append(ms1MatchResult,ppm)
  names(ms1MatchResult)[length(ms1MatchResult)] <- "ppm"
  #加入ppm属性
  object@ms1MatchResult <- ms1MatchResult
  #存储结果
  return(object)
}

#' @title Match MS1 with the CL-database.
#' @param object a findmyCL object.
#' @param ppm numeric(1),default is 5.
#' @param database_name character(1)
#' @seealso [findmyCL::matchDatabase()]
#' @return dataframe(1)
#' @export
matchMS1_main <- function(object , ppm = 5 , database_name = "CL"){
  data <-  object@xcms@peaks
  #需要计算的数据
  data_deal <- t(data)
  vector <- lapply(seq_len(ncol(data_deal)), function(i) data_deal[,i])
  #将数据框转换为列表便于多线程
  paste0("Matching with " ,database_name ," ...") %>%
    message()
  pb <- dplyr::progress_estimated(length(vector))
  #设置进度条
  options (warn = -1)
  #显示提示信息
  result <- purrr::map(.x = vector,.f = findmyCL::matchDatabase, vector_column = 1 ,ppm = ppm,database = findmyCL:::MS1_database[[database_name]],database_column = 1,pb = pb)
  #批量计算结果,记得这里参数要传入pb
  cat("\n")
  #换行
  result <- findmyCL::deleteNULL(list = result)
  #去除空值
  result <- plyr::rbind.fill(result)
  #计算配对结果

  colnames_length <- length(colnames(data))
  #计算原本数据列名长度
  result_name_length <- length(colnames(result))
  #计算结果的列名长度
  colnames(result)[(colnames_length + 1):result_name_length] <- c("aditive form","Chain Length:Δ","Formula","Oxform")
  colnames(result)[1:colnames_length] <- colnames(data)
  #改结果的列名
  options (warn = 1)
  #显示提示信息
  return(result)
}
