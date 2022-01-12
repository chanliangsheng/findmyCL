#' @title Delete NA in a vector
#' @description Delete NA in a vector.
#' @param vector vector(1)
#' @return vector(1)
#' @export
deleteNA <- function(vector){
  if (length(which(is.na(vector))) == 0) {
    return(vector)
  }
  #如果没有NA，则返回原本的向量
  if (length(which(is.na(vector))) > 0) {
    vector <- vector[-which(is.na(vector))]
    return(vector)
  }
  #如果有NA，则删除NA
}
