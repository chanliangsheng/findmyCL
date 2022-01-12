#' @title Turn a dataframe into list by row
#' @param dataframe dataframe(1)
#' @return list(1)
#' @export
turnDataframeList <- function(dataframe){
  data_deal <- t(dataframe)
  list <- lapply(seq_len(ncol(data_deal)), function(i) data_deal[,i]) %>%
    return()
  #转换为列表并返回
}
