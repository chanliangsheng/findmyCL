#' @title Delete the NULL element in a list
#' @description Delete the NULL element in a list.
#' @param list list(1)
#' @export
deleteNULL <- function(list){
  if (length(which(sapply(list, is.null))) != 0) {
    list <-  list[-which(sapply(list, is.null))]
    #去除NULL
  }
  #如果有NULL值，则将其去除，如果没有NULL值，则不对list进行操作
  return(list)
}
