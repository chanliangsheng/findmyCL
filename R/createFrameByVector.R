#' @title Creating a dataframe by vector
#' @description Decompose a vector.Every n elements in this vector constitute a row in the new dataframe.
#' @param vector vector(1)
#' @param n numeric(1)
#'
#' @return dataframe(1)
#' @export
createFrameByVector <- function(vector , n){
  result <- c()
  for (i in 1:n) {
    result <- rbind(result,vector)
  }
  result <- as.data.frame(result)
  #转换为数据框
  row.names(result) <- 1:n
  return(result)
}
