#' @title Checking which MS1 have MS2,from[matchMS1()]
#' @description Checking if a MS1 have corresponding MS2,use this function after [matchMS1()].A MS1 may have more than 1 MS2.The result is stored in ms2CheckResult.
#' @param object a findmyCL object
#' @param ppm numeric(1) defining the maximal tolerated m/z deviation in consecutive scans in parts per million (ppm) for the initial ROI definition when matching with the precursorMz of MS2,default is 5.
#' @return a findmyCL object
#' @export
checkMS2 <- function(object , ppm = 5){
  dealing_name <- as.list(names(object@ms1MatchResult))
  dealing_name <- head(dealing_name,-1)
  #转换为列表方便迭代操作,去除ppm属性
  ms2CheckResult <- purrr::map(.x = dealing_name,.f = findmyCL::checkMS2_main,object = object,ppm = ppm)
  #迭代运算不同的一级配对结果
  noMS2 <- purrr::map(.x = dealing_name,.f = findmyCL::checknoMS2_main,object = object,ppm = ppm)
  #迭代运算不同的一级配对结果
  names(ms2CheckResult) <- dealing_name
  names(noMS2) <- dealing_name
  #重命名
  ms2CheckResult <- append(ms2CheckResult,ppm)
  names(ms2CheckResult)[length(ms2CheckResult)] <- "ppm"
  noMS2 <- append(noMS2,ppm)
  names(noMS2)[length(noMS2)] <- "ppm"
  #加入ppm属性
  object@ms2CheckResult <- ms2CheckResult
  object@noMS2 <- noMS2
  #加入结果
  return(object)
}
#识别一个对象中的matchMS1的结果是否有MS2的配对结果

#' @title Check every MS1 with its MS2 in a MS1-result
#' @param object a findmyCL object
#' @param ppm numeric(1) defining the maximal tolerated m/z deviation in consecutive scans in parts per million (ppm) for the initial ROI definition when matching with the precursorMz of MS2,default is 5.
#' @param ms1matchresult_name character(1)
#' @import MSnbase
#' @return list(1)
#' @export
checkMS2_main <- function(object , ppm = 5 , ms1matchresult_name = "CL"){
  nrow <- length(object@MS2@assayData)
  ncol <- 3
  database <- as.data.frame(matrix(nrow = nrow,ncol = ncol))
  names(database) <- c("name","precursorMz","rt")
  #创建数据库框架
  MS2 <- object@MS2@assayData %>%
    as.list()
  database$name <- names(MS2)
  precursorMz <- MS2 %>%
    purrr::map(function(x) x@precursorMz) %>%
    as.numeric()
  rt <- MS2 %>%
    purrr::map(function(x) x@rt) %>%
    as.numeric()
  database$precursorMz <- precursorMz
  database$rt <- rt
  #创建MS2

  data <-  object@ms1MatchResult[[ms1matchresult_name]]
  #需要计算的数据
  data_deal <- t(data)
  vector <- lapply(seq_len(ncol(data_deal)), function(i) data_deal[,i])
  #将数据框转换为列表便于多线程
  paste0("Checking with " ,ms1matchresult_name ," ...") %>%
    message()
  pb <- dplyr::progress_estimated(length(vector))
  #设置进度条
  result <- purrr::map(.x = vector,.f = findmyCL::searchMS2,object = object , vector_mz_column = 1 , vector_rtmin_column = 5,vector_rtmax_column = 6, ppm = ppm , MS2 = database , MS2_column = 2 , pb = pb)
  #批量处理
  cat("\n")
  #换行
  result <- findmyCL::deleteNULL(result)
  #去除空值
  names(result) <- 1:length(result)
  #重命名结果
  message("Done!")
  return(result)
}
#为每一个matchMS1的结果比如“CL”中的一级峰识别是否存在二级，并整合为列表

#' @title Check which MS1 has no MS2 in every MS1-match-result
#' @param object a findmyCL object
#' @param ppm numeric(1) defining the maximal tolerated m/z deviation in consecutive scans in parts per million (ppm) for the initial ROI definition when matching with the precursorMz of MS2,default is 5.
#' @param ms1matchresult_name character(1)
#' @return list(1)
#' @export
checknoMS2_main <- function(object , ppm = 5 , ms1matchresult_name = "CL"){
  nrow <- length(object@MS2@assayData)
  ncol <- 3
  database <- as.data.frame(matrix(nrow = nrow,ncol = ncol))
  names(database) <- c("name","precursorMz","rt")
  #创建数据库框架
  MS2 <- object@MS2@assayData %>%
    as.list()
  database$name <- names(MS2)
  precursorMz <- MS2 %>%
    purrr::map(function(x) x@precursorMz) %>%
    as.numeric()
  rt <- MS2 %>%
    purrr::map(function(x) x@rt) %>%
    as.numeric()
  database$precursorMz <- precursorMz
  database$rt <- rt
  #创建MS2

  data <-  object@ms1MatchResult[[ms1matchresult_name]]
  #需要计算的数据
  data_deal <- t(data)
  vector <- lapply(seq_len(ncol(data_deal)), function(i) data_deal[,i])
  #将数据框转换为列表便于多线程
  paste0("Checking no MS2 with " ,ms1matchresult_name ," ...") %>%
    message()
  pb <- dplyr::progress_estimated(length(vector))
  #设置进度条
  result <- purrr::map(.x = vector,.f = findmyCL::searchNoMS2,object = object , vector_mz_column = 1 , vector_rtmin_column = 5,vector_rtmax_column = 6, ppm = ppm , MS2 = database , MS2_column = 2 , pb = pb)
  #批量处理
  cat("\n")
  #换行
  result <- findmyCL::deleteNULL(result)
  #去除空值
  names(result) <- 1:length(result)
  #重命名结果
  message("Done!")
  return(result)

}
#为每一个matchMS1的结果比如“CL”中不存在二级的一级，并整合为列表

#' @title Search which MS1 have its MS2
#' @param object a findmyCL object
#' @param vector vector(1)
#' @param vector_mz_column numeric(1)
#' @param vector_rtmin_column numeric(1)
#' @param vector_rtmax_column numeric(1)
#' @param ppm numeric(1),defining the maximal tolerated m/z deviation in consecutive scans in parts per million (ppm) for the initial ROI definition when matching with the precursorMz of MS2,default is 5.
#' @param MS2 dataframe(1)
#' @param MS2_column numeric(1)
#' @param pb process bar
#' @return list(1)
#' @export
searchMS2 <- function(object , vector , vector_mz_column , vector_rtmin_column,vector_rtmax_column , ppm = 5 , MS2 , MS2_column , pb){
  pb$tick()$print()
  #显示进度条
  min_match <- MS2[,MS2_column] - MS2[,MS2_column]*ppm/1000000
  max_match <- MS2[,MS2_column] + MS2[,MS2_column]*ppm/1000000
  #计算数据库配对的最小值和最大值
  match_number_mz <- which((as.numeric(vector[vector_mz_column]) > min_match) & (as.numeric(vector[vector_mz_column]) < max_match))
  #满足一级mz在二级的前体离子的波动范围的二级
  match_number_rt <- which((as.numeric(vector[vector_rtmin_column]) < MS2$rt) &(as.numeric(vector[vector_rtmax_column]) > MS2$rt))
  #满足二级的保留时间在一级保留时间的波动范围内的二级
  match_number <- intersect(match_number_mz,match_number_rt)
  #求交集，配对
  if (length(match_number) == 0) {
    return(NULL)
  }
  #如果没有二级，返回空值
  if (length(match_number) > 0) {
    MS2_index <- MS2[match_number,][1]
    MS2_index <- as.list(t(MS2_index))
    #获取对应的MS2的索引
    result <- purrr::map(.x = MS2_index,.f = findmyCL::catchMS2,object = object)
    #提取匹配到的二级的结果

    store <- list()
    vector <- findmyCL::deleteNA(vector = vector)
    #将一级配对结果中多余的NA去除
    MS1 <- findmyCL::turnVectorDataframe(vector = as.data.frame(t(as.matrix(vector))))
    #把向量转换为数据框,输入的需要是一个数据框，所以需要转变为一行的数据框
    store <- findmyCL::fappend(list1 = store,list2 = MS1)
    store <- append(x = store,values = result)
    #整合到store中，第一个元素是一级配对的结果，剩下的都是对于的二级

    names(store) <- paste0("MS2_",0:(length(store) - 1))
    names(store)[1] <- "MS1"
    #重命名结果
    return(store)
  }
}
#识别单个一级峰是否存在二级，并整合为列表，可能会有多个二级


#' @title Search which MS1 have no MS2
#' @param object a findmyCL object
#' @param vector vector(1)
#' @param vector_mz_column numeric(1)
#' @param vector_rtmin_column numeric(1)
#' @param vector_rtmax_column numeric(1)
#' @param ppm numeric(1),defining the maximal tolerated m/z deviation in consecutive scans in parts per million (ppm) for the initial ROI definition when matching with the precursorMz of MS2,default is 5.
#' @param MS2 dataframe(1)
#' @param MS2_column numeric(1)
#' @param pb process bar
#' @return list(1)
#' @export
searchNoMS2 <- function(object , vector , vector_mz_column , vector_rtmin_column,vector_rtmax_column , ppm = 5 , MS2 , MS2_column , pb){
  pb$tick()$print()
  #显示进度条
  min_match <- MS2[,MS2_column] - MS2[,MS2_column]*ppm/1000000
  max_match <- MS2[,MS2_column] + MS2[,MS2_column]*ppm/1000000
  #计算数据库配对的最小值和最大值
  match_number_mz <- which((as.numeric(vector[vector_mz_column]) > min_match) & (as.numeric(vector[vector_mz_column]) < max_match))
  #满足一级mz在二级的前体离子的波动范围的二级
  match_number_rt <- which((as.numeric(vector[vector_rtmin_column]) < MS2$rt) &(as.numeric(vector[vector_rtmax_column]) > MS2$rt))
  #满足二级的保留时间在一级保留时间的波动范围内的二级
  match_number <- intersect(match_number_mz,match_number_rt)
  #求交集，配对
  if (length(match_number) > 0) {
    return(NULL)
  }
  #如果有二级，则返回NULL
  if (length(match_number) == 0) {
    vector <- findmyCL::deleteNA(vector = vector)
    #将一级配对结果中多余的NA去除
    return(vector)
  }
}

#' @title Catch MS2
#' @param object a findmyCL object
#' @param MS2_name character(1)
#' @return Spectrum2
#' @export
catchMS2 <- function(object , MS2_name){
  object@MS2@assayData[[MS2_name]] %>%
    return()
}
