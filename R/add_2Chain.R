#' @title Adding two Chain Length:Δ
#' @description This function is used to add 2 chain:Δ(it can be 2 vector).Like add "10:1" + "11:1" = "21:2".
#' @param vector1 vector(1),like "58:1"
#' @param vector2 vector(1),like "59:2"
#' @return vector(1)
#' @export
add_2Chain <- function(vector1 , vector2){
  chain_sum <- as.numeric(as.data.frame(strsplit(vector1,":"))[1,]) + as.numeric(as.data.frame(strsplit(vector2,":"))[1,])
  #求chain相加的结果
  double_sum <- as.numeric(as.data.frame(strsplit(vector1,":"))[2,]) + as.numeric(as.data.frame(strsplit(vector2,":"))[2,])
  #求双键相加的结果

  result <- paste0(as.character(chain_sum),":",as.character(double_sum)) %>%
    return()
  #返回结果
}
