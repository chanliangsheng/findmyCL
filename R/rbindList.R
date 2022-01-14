#' @title Bind all dataframe in a list with same colnames.
#' @param list list(1)
#' @return dataframe(1)
#' @export
rbindList <- function(list){
  colname <- colnames(list[[1]])
  #将原本的列名存储起来
  result <- purrr::map(.x = list , .f = findmyCL::rbindList_main) %>%
    plyr::rbind.fill()
  #改变所有数据框的列名,合并所有数据框
  colnames(result) <- colname
  #将存储的列名返回
  return(result)
  #返回结果
}

#' @title Change every colnames of dataframe in the list.
#' @param dataframe dataframe(1)
#' @return dataframe(1)
#' @export
rbindList_main <- function(dataframe){
  colnames(dataframe) <- 1:length(colnames(dataframe))
  #重命名列名方便行合并
  return(dataframe)
  #返回结果
}
