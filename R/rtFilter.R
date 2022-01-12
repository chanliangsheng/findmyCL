#' @title Split every MS2matchresult into half
#' @description This is a function to split MS2matchresult into half.Half have only one cardiolipin,the other have more than tow match cardiolipin.The former is used to select which is the right cardiolipin in the latter half.
#' @param object a findmyCL class
#' @return a findmyCL class
#' @export
rtfilter <- function(object){
  deal_names <- head(names(object@ms2MatchResult) , -1)
  #去除最后一个ppm的属性的名字
  result <- purrr::map(.x = as.list(deal_names) , .f = findmyCL::rtfilter_main , object = object)
  #迭代计算
  names(result) <- deal_names
  #重命名结果
  result <- append(x =  result , values =  tail(object@ms2MatchResult , 1))
  #追加ppm属性
  return(result)
}

#' @title Split a MS2match&splice result into half
#' @description Delete the redundant Cardiolipin in MS1(a MS1 may have more than 1 possible Cardiolipin and it may can not be spliced into Cardiolipin).Delete redundant splicePA result in every MS2(some splicePA result may can not be spliced by its FA in its MS2).Work for [findmyCL::rtfilter()]
#' @param object a findmyCL class
#' @param deal_name character(1) the name of MS2matchresult,like "CL"
#' @return list(1)
#' @export
rtfilter_main <- function(object , deal_name = "CL"){
  delete_result <- purrr::map(.x = object@ms2MatchResult[[deal_name]] ,.f = delete_no_FA_MS1)
  #去除MS1中不能由PA或者FA拼接成心磷脂的可能心磷脂

  rtseed <- purrr::map(.x = delete_result , .f = findmyCL::choose_list , mode = 1) %>%
    findmyCL::deleteNULL()
  #寻找1个一级对应1个心磷脂的bigMS1
  rtscreen <- purrr::map(.x = delete_result , .f = findmyCL::choose_list , mode = 2) %>%
    findmyCL::deleteNULL()
  #寻找1个一级对应多个心磷脂的bigMS1
  if (length(rtseed) > 0) {
    names(rtseed) <- 1:length(rtseed)
  }
  #重命名结果
  if (length(rtscreen) > 0) {
    names(rtscreen) <- 1:length(rtscreen)
  }
  #重命名结果

  result_list <- list(rtseed , rtscreen)
  #整合结果
  names(result_list) <- c("rtseed","rtscreen")
  #重命名结果

  return(result_list)
}
#将二级拼接结果中一级对应1个和多个的心磷脂分开，因为1个一级对应多个心磷脂的话要筛选

#' @title Choose 1-MS1 vs 1 Cardiolipin and 1-MS1 vs more Cardiolipin.
#' @param list list(1),a bigMS1 contain MS1,MS2,splice-Cardiolipin-result.
#' @param mode numeric(1),modes used to select 1 vs 1 or 1 vs more than 1.
#' @return list(1)
#' @export
choose_list <- function(list , mode = 1){
  if (mode == 1) {
    if (length(list$MS1[,1]) == 1) {
      return(list)
    }
  }
  #mode1：MS1的结果是1个一级对应1个心磷脂

  if (mode == 2) {
    if (length(list$MS1[,1]) > 1) {
      return(list)
    }
  }
  #mode2：MS1的结果是1个一级对应多个心磷脂
}
#选择有对应多个心磷脂和单个心磷脂的一级

#' @title Delete redundant Cardiolipin.
#' @description Delete redundant Cardiolipin and splicePA-result in MS1 which can not spliced by PA&FA.
#' @param bigMS1 list(1),a bigMS1 contain MS1,MS2,splice-Cardiolipin-result.
#' @return list(1),a bigMS1 contain MS1,MS2,splice-Cardiolipin-result.
#' @export
delete_no_FA_MS1 <- function(bigMS1){
  if (length(bigMS1$MS1[,1]) == 1) {
    return(bigMS1)
  }
  #如果MS1中只有一种配对心磷脂的可能，则返回该bigMS1
  deal_list <- bigMS1[-1]
  #去除MS1,提取所有的MS2
  all_chain <- purrr::map(.x = deal_list , .f = findmyCL::get_all_chain) %>%
    unlist() %>%
    unique()
  #提取所有MS2中的配对成功的所有Chain Length:Δ
  delete_redundant_PA_result <- purrr::map(.x = deal_list , .f = findmyCL::delete_redundant_PA , chain = all_chain)
  #去除二级拼接结果中多余的splicePA（例如：即一级有对应的两个心磷脂，但是只有一个能由PA和FA拼接，另一个只能找到PA，找不到FA，要将这个找不到FA的PA的拼接结果去除）
  bigMS1 <- append(x = bigMS1[1] , delete_redundant_PA_result)
  #将bigMS1中的一级的结果和去除多余拼接PA的结果整合为新的bigMS1

  MS1matchresult <- match(all_chain , bigMS1$MS1$`Chain Length:Δ`) %>%
    sort()
  #匹配哪些心磷脂能够被拼接出来

  bigMS1$MS1 <- bigMS1$MS1[MS1matchresult , ]
  #将可能的心磷脂赋值回原来的位置


  return(bigMS1)
}
#检验一级配对的结果里面的哪些心磷脂可以被拼接出来，返回一个bigMS1

#' @title Get all Chain Length:Δ in a list.
#' @description Extract a Chain Length:Δ in a splice-result(a list).
#' @param list list(1),contain splice-result and Chain Length:Δ.
#' @return character(1)
#' @seealso [findmyCL::get_all_chain()]
#' @export
get_chain <- function(list){
  list$`Chain Length:Δ` %>%
    return()
}
#提取一个列表中的所有Chain Length:Δ

#' @title Get all Chain Length:Δ in a MS2.
#' @description Extract all Chain Length:Δ from splicePA-result && spliceFA-result in a MS2.
#' @param MS2 list(1)
#' @return matrix(1)
#' @export
get_all_chain <- function(MS2){
  if (is.null(MS2$spliceFA)) {
    chain <- purrr::map(.x = MS2[["splicePA"]] , .f = findmyCL::get_chain) %>%
      as.matrix() %>%
      return()
  }
  #如果没有拼接成FA，只拼接成PA,则提取所有splicePA中的Chain Length:Δ

  if (is.null(MS2$spliceFA) == FALSE) {
    chain <- purrr::map(.x = MS2[["spliceFA"]] , .f = findmyCL::get_chain) %>%
      as.matrix() %>%
      return()
  }
  #如果拼接成了FA，只拼接成PA,则提取所有spliceFA中的Chain Length:Δ
}
#提取一个MS2中FA配对成功的所有Chain Length:Δ

#' @title Delete redundant splicePA-result.
#' @description Use Chain Length:Δ in spliceFA-result(it exist mean a Cardiolipin can be spliced by PA which also can be spliced by FA).But some splicePA-result can not be spliced by FA(it has not be deleted yet).So use this function to delete these splicePA-result.Work for [findmyCL::delete_no_FA_MS1()]
#' @param MS2 list(1)
#' @param chain vector(1) or character(1)
#' @return list(1)
#' @export
delete_redundant_PA <- function(MS2 , chain){
  if (class(MS2[["PA"]]) == "numeric") {
    return(MS2)
  }
  #如果MS2中的PA的类型是数值型而不是数据框，说明PA中就是0，即没有找到PA
  #如果MS2寻找PA和FA的过程中没有找到PA，则不删除多余的PA（因为根本没有PA拼接的结果）
  splicePA <- purrr::map(.x = MS2[["splicePA"]] , .f = findmyCL::delete_redundant_PA_main , chain = chain) %>%
    findmyCL::deleteNULL()
  #去除多余的splicePA中的结果
  names(splicePA) <- 1:length(splicePA)
  names(MS2$spliceFA) <- 1:length(MS2$spliceFA)
  #重命名结果
  MS2$splicePA <- splicePA
  #将新的结果赋值在原来结果上
  return(MS2)
  #返回该二级
}
#这里的二级指的是包含拼接结果的二级！！！！

#' @title Delete redundant splicePA-result.
#' @description Work for [findmyCL::delete_redundant_PA()]
#' @param splicePA_result list(1),contain splicePA,Chain Length:Δ,oxygen
#' @param chain vector(1)
#' @return list(1)
#' @export
delete_redundant_PA_main <- function(splicePA_result  , chain){
  match <- length(which(splicePA_result[["Chain Length:Δ"]] == chain))
  #判断这个splicePA的结果的Chain Length:Δ与spliceFA的结果(即chain，包含了FA中的所有Chain Length:Δ)是否对应
  if (match == 0) {
    splicePA_result <- NULL
    return(splicePA_result)
  }
  #如果与FA的拼接结果的Chain Length:Δ不符合，则返回NULL
  if (match != 0) {
    return(splicePA_result)
  }
  #如果与FA的拼接结果的Chain Length:Δ符合，则返回原本的splicePA_result
}
#对splicePA中的每个结果进行检查
