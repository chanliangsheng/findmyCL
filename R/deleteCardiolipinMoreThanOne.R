#' @title Delete redundant Cardiolipin
#' @description Some MS1 have more than 1 corresponding Cardiolipin,but it should only one  Cardiolipin in a MS1,this function use seed which is actually a Cardiolipin get only one Cardiolipin named rtseed to choose which Cardiolipin is the right Cardiolipin in MS1 which get more than 2 Cardiolipin.Use retention time rule to filter some chain:Δ.The longer the chain,the longer the retention time.The greater the unsaturation,The shorter the retention time.
#' @param object a findmyCL object
#' @return a findmyCL object
#' @export
deleteCardiolipinMoreThanOne <- function(object){

}

#' @title Delete redundant Cardiolipin for a single MS2matchresult
#' @param object a findmyCL object
#' @param ms2matchresult_name character(1)
#' @return list(1)
#' @export
deleteCardiolipinMoreThanOne_main <- function(object , ms2matchresult_name = "CL"){
  result <- purrr::map(.x = object@ms2MatchResult[[ms2matchresult_name]][["rtseed"]] , .f = findmyCL::getAllChainRt) %>%
    rbind.fill()
  #提取rt种子中的所有rt,chain,Δ放在一个数据框中

}

#' @title Get all rt,chain,Δ in a bigMS1.
#' @description Get all rt,chain,Δ in a bigMS1.A bigMS1 is consist of MS1 and some MS2.Take rt,chain,Δ
#' @param bigMS1 list(1)
#' @return dataframe(1)
#' @export
getAllChainRt <- function(bigMS1){
  store <- matrix(nrow = 1 , ncol = 3) %>%
    as.data.frame()
  colnames(store) <- c("rt","chain","Δ")
  #重命名
  chain_Unsaturation <- strsplit(bigMS1$MS1$`Chain Length:Δ` , ":")[[1]] %>%
    as.vector()
  #提取chain 和 Δ
  store$rt <- as.numeric(bigMS1$MS1$rt)
  store$chain <- as.numeric(chain_Unsaturation[1])
  store$Δ <- as.numeric(chain_Unsaturation[2])
  #提取rt、chain、Δ
  return(store)
  #返回结果
}
#bigMS1是包含MS1和多个MS2的集合
