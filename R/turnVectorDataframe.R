#' @title Turn a vector to a dataframe
#' @description A MS1 may has only 1 corresponding Cardiolipin.But some MS1 may has more than 1 corresponding Cardiolipins.So this function is used to turn these 1 vs much Cardiolipin(vector) into a dataframe.
#' @param vector vector(1)
#' @return dataframe(1)
#' @export
turnVectorDataframe <- function(vector){
  if (length(vector) == 15) {
    return(vector)
  }
  #如果只有一个匹配心磷脂的结果，则直接返回这个向量
  vector_length <- length(vector)
  #求向量长度
  factor <- gl((vector_length - 11)/4 , 4)
  split_result <- split(x = as.vector(as.matrix(vector))[12:vector_length] , f = factor)
  #拼接结果
  split_result <- as.data.frame(t(as.data.frame(split_result)))
  #拼接在一起
  options(warn = -1)
  #不显示warning
  split_result <- cbind(vector[1,1:11],split_result)
  #拼接vector和它的心磷脂
  options(warn = 1)
  #显示warning
  colnames(split_result) <- colnames(vector)[1:15]
  #重命名
  return(split_result)
}
