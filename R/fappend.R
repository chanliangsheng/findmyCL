#' @title Append all elements in list2 to list1.
#' @description A improvement of [append()].If you use test <- append(x = list1 , values = list2).The result is the list2 is added to the list1.But what we want is add whole elements in list2 to list1 rather than adding a list2 to list1.So this function
#' @param list1 the list appended
#' @param list2 the list you choose to append
#' @return list(1)
#' @export
fappend <- function(list1,list2){
  list_help <- list(list(1))
  #新建一个列表用来帮助转移
  list_help[[1]] <- list2
  #将list2放于新建的这个列表的第一个元素中
  names(list_help) <- deparse(substitute(list2))
  #改变list2的名称为原来的名称
  return(append(list1,list_help))
  #返回将list2append到list1的结果
}
