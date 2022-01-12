#' @title Matching dataframe with a database
#' @description Matching dataframe with a database consist of list.This function is used in [findmyCL::matchMS1()].
#' @param vector vector(1)
#' @param vector_column numeric(1) which column you want to compare
#' @param ppm The floating value you can tolerate most
#' @param database dataframe(1)
#' @param database_column numeric(1) which column in database you choose to compare
#' @return a dataframe
#' @export
matchDatabase <- function(vector , vector_column , ppm = 5 , database , database_column ,pb){
  min_match <- database[,database_column] - database[,database_column]*ppm/1000000
  max_match <- database[,database_column] + database[,database_column]*ppm/1000000
  #计算数据库配对的最小值和最大值
  match_number <- which((vector[vector_column] > min_match) & (vector[vector_column] < max_match))
  #计算满足数据库中的哪些值
  pb$tick()$print()
  #显示进度条
  if (length(match_number) == 0) {
    return(NULL)
  }
  #如果没有配对成功，则返回空值
  if (length(match_number) > 0){
    database_match <-  database[match_number,-database_column]
    #取数据库中除了用来配对的数据，其他的用来合并到结果中
    database_match <- as.data.frame(t(as.matrix(as.vector(t(as.matrix(database_match))))))
    #转变为一行的数据框
    vector <- as.data.frame(t(as.matrix(vector)))
    #转变为数据框才可以列合并
    match_result <- cbind(vector,database_match)
    #合并结果
    names(match_result) <- 1:length(match_result)
    #改变列名方便合并
    return(match_result)
  }
}
