#' @title Splice PA and FA from a MS2 and MS1 in DLCL
#' @description DLCL:Double hydrolyzed cardiolipin(have 1 PA).A function to search PA and FA from MS2.If this MS1 can not search its 2 PA,delete this MS1(this MS1 may not have just one corresponding Cardiolipin).If this MS2 can search its 2 PA but can not search corresponding 2 FA or it can just match one FA,delete this MS2.
#' @import purrr
#' @import dplyr
#' @param MS2 list(1),including MS2 information,PA&FA search result
#' @param MS1 dataframe(1),including mz,rt,rtmin,rtmax ...
#' @return list(1),a new MS2,adding PA&FA splicing result
#' @export
splicePAFA_DLCL <- function(MS2 , MS1){
  if (is.null(MS2)) {
    return(NULL)
  }
  #如果传入的MS2是空值，则直接返回空值，因为这代表着传入的MS2即没有PA也没有FA
  havePA <- length(MS2$PA)
  haveFA <- length(MS2$FA)
  #计算MS2是没有PA还是没有FA，还是两者都有
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if ((havePA == 1) & (haveFA == 5)){
    splicePA_list <- purrr::pmap(.l = list(as.list(MS1$`Chain Length:Δ`) , as.list(MS1$Oxform)) , .f = findmyCL::splice2FA_DLCL , FA = MS2$FA) %>%
      findmyCL::deleteNULL()
    #进行PA的拼接，去除NULL
    if (length(splicePA_list) == 0) {
      return(NULL)
    }
    #如果这个二级的FA都无法拼接成DLCL，则将这个二级去除（即返回NULL）
    MS2 <- findmyCL::fappend(list1 = MS2 , list2 = spliceFA_list)
    #追加到MS2中
    names(MS2) <- c("MS2","PA","FA","spliceFA")
    #重命名追加结果
    return(MS2)
  }
  #如果没有PA，只有FA，则只进行FA的拼接，2个FA能组成DLCL就好
  if ((haveFA == 1) & (havePA == 5)){
    splicePA_list <- purrr::pmap(.l = list(as.list(MS1$`Chain Length:Δ`) , as.list(MS1$Oxform)) , .f = findmyCL::splicePA_DLCL , PA = MS2$PA) %>%
      findmyCL::deleteNULL()
    #进行PA的拼接，去除NUL
    if (length(splicePA_list) == 0) {
      return(NULL)
    }
    #如果这个二级的PA都无法拼接成DLCL，则将这个二级去除（即返回NULL）
    names(splicePA_list) <- 1:length(splicePA_list)
    #重命名拼接PA的结果
    MS2 <- findmyCL::fappend(list1 = MS2 , list2 = splicePA_list)
    #追加到MS2中
    names(MS2) <- c("MS2","PA","FA","splicePA")
    #重命名追加结果
    return(MS2)
    #返回结果
  }
  #假如没有FA，则只拼接PA
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if ((haveFA == 5) & (havePA == 5)){
    splicePA_list <- purrr::pmap(.l = list(as.list(MS1$`Chain Length:Δ`) , as.list(MS1$Oxform)) , .f = findmyCL::splicePA_DLCL , PA = MS2$PA) %>%
      findmyCL::deleteNULL()
    #进行PA的拼接，去除NULL,splicePA_list包括多个一级对应的拼接结果，是一个列表，每个属性包括每个一级对应的结果
    if (length(splicePA_list) == 0) {
      return(NULL)
    }
    #如果这个二级的PA都无法拼接成DLCL，则将这个二级去除（即返回NULL）
    names(splicePA_list) <- 1:length(splicePA_list)
    #重命名拼接PA的结果
    spliceFA_list <- purrr::map(.x = splicePA_list , .f = findmyCL::splicePA2FA_DLCL , FA = MS2$FA) %>%
      findmyCL::deleteNULL()
    #用PA的拼接结果拼接2个FA，并且找另外一个FA能一起拼接成DLCL
    if (length(spliceFA_list) == 0) {
      return(NULL)
    }
    #如果不能拼接成DLCL，则返回NULL
    MS2 <- findmyCL::fappend(list1 = MS2 , list2 = splicePA_list)
    MS2 <- findmyCL::fappend(list1 = MS2 , list2 = spliceFA_list)
    #追加到MS2中
    names(MS2) <- c("MS2","PA","FA","splicePA","spliceFA")
    #重命名追加结果
    return(MS2)
  }
  #如果既有PA，又有FA，则先进行PA的拼接，再进行FA的拼接
}


#' @title Search which PA can be a DLCL.
#' @description A PA can be a DLCL if its Chain Length:Δ and oxform can match this DLCL.
#' @param PA dataframe(1)
#' @param chain_double character(1)
#' @param oxygen numeric(1)
#' @return list(1),including splicePA-result,its corresponding Chain Length:Δ && oxygen
#' @export
splicePA_DLCL <- function(PA , chain_double , oxygen){
  chain_match <- which(PA$`Chain Length:Δ` == chain_double)
  oxygen_match <- which(PA$Oxform == oxygen)
  #求哪些PA的Chain Length:Δ与一级的Chain Length:Δ相同；oxygen与一级的oxygen相同,返回行号
  match_result <- intersect(chain_match , oxygen_match)
  #求交集，Chain Length:Δ和oxygen都与一级的相同的行号
  if (length(match_result) == 0) {
    return(NULL)
  }
  #如果没有配对成功，返回NULL
  PA_match_result <- PA[match_result,]
  #配对成功的取出
  rownames(PA_match_result) <- 1:length(PA_match_result[,1])
  #重命名行号
  result_list <- list(PA_match_result , chain_double , oxygen)
  #生成列表为结果
  names(result_list) <- c("PA","Chain Length:Δ","oxygen")
  #重命名结果
  return(result_list)
  #返回结果
}

#' @title Splice 2FA into a PA.
#' @description In DLCL,we find a PA first.Then,search which FA can compose this PA.
#' @param splicePA list(1),including PA,Chain Length:Δ,oxygen.
#' @param FA dataframe(1)
#' @return list(1)
#' @export
splicePA2FA_DLCL <- function(splicePA , FA){
  list <- findmyCL::turnDataframeList(dataframe = splicePA[["PA"]])
  #将数据框转换为列表便于多线程
  spliceFA_result <- purrr::map(.x = list , .f = findmyCL::splice2FA_MLCL_main , FA = FA)
  #求由FA拼接成PA的结果
  spliceFA_result <- findmyCL::deleteNULL(spliceFA_result)
  #去除NULL
  if (length(spliceFA_result) == 0) {
    return(NULL)
  }
  #如果FA没有配对成功，则返回空值
  spliceFA_result <- findmyCL::rbindList(list = spliceFA_result)
  #从列表转换为数据框，并且将所有结果行整合
  store <- list(spliceFA_result , splicePA[["Chain Length:Δ"]] , splicePA[["oxygen"]])
  #存储FA拼接结果。。。
  names(store) <- c("FA" , "Chain Length:Δ" , "oxygen")
  #重命名结果
  return(store)
}

#' @title Splice 2FA into a DLCL
#' @description Splice 2FA into a DLCL which have no PA.
#' @param FA dataframe(1)
#' @param chain_double character(1)，like "14:0".
#' @param oxygen numeric(1)
#' @return list(1)
#' @export
splice2FA_DLCL <- function(FA , chain_double , oxygen){
  FA_num <- length(FA[,1])
  #求一共有多少个FA
  FA_combination <- expand.grid(rep(list(1 : FA_num), 2)) %>%
    apply(1 , sort , decreasing = T) %>%
    t() %>%
    as.data.frame() %>%
    dplyr::distinct()
  #计算所有的组合可能(行号的组合)
  chain_add_result <- findmyCL::add_2Chain(vector1 = FA[FA_combination[ , 1] , 4] , FA[FA_combination[ , 2] , 4])
  #将所有组合可能的chain:double相加，2行相加(FA中的chain:double相加)
  oxygen_add_result <- FA[FA_combination[ , 1] , 5] + FA[FA_combination[ , 2] , 5]
  #将所有组合可能的氧原子个数相加，2行相加(FA中的oxygen相加)
  matchresult <- intersect(which(chain_add_result == chain_double) , which(oxygen_add_result == oxygen))
  #取即与chain：double相等又与氧原子个数相等的行号
  if (length(matchresult) == 0) {
    return(NULL)
  }
  #如果配对没有成功，则返回空值
  FA_combination_result <- FA_combination[matchresult , ]
  #将配对成功结果的行号保存
  catch_result <- purrr::pmap(.l = list(row1 = list(FA_combination_result[,1]) , row2 = list(FA_combination_result[,2])), .f = findmyCL::catchDataframe2Col , dataframe = FA)
  #将所有符合的行提取
  colnames(catch_result[[1]]) <- rep(colnames(FA),2)
  #重命名FA拼接结果的行名
  catch_result <- append(x = catch_result , values = chain_double)
  catch_result <- append(x = catch_result , values = oxygen)
  names(catch_result) <- c("FA","Chain Length:Δ" , "oxygen")
  #追加属性，重命名结果
  return(catch_result)
  #返回结果

}
