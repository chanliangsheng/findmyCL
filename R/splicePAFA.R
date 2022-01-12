#' @title Splice PA and FA from a MS2 and MS1 in CL
#' @description A function to search PA and FA from MS2.If this MS1 can not search its 2 PA,delete this MS1(this MS1 may not have just one corresponding Cardiolipin).If this MS2 can search its PA but can not searth corresponding 2 FA or it can just match one FA,delete this MS2.
#' @import gtools
#' @import purrr
#' @import furrr
#' @import dplyr
#' @importFrom tidyr crossing
#' @param MS2 list(1),including MS2 information,PA&FA search result
#' @param MS1 dataframe(1),including mz,rt,rtmin,rtmax ...
#' @return list(1),a new MS2,adding PA&FA splicing result
#' @export
splicePAFA <- function(MS2 , MS1){
  if (is.null(MS2)) {
    return(NULL)
  }
  #如果传入的MS2是空值，则直接返回空值，因为这代表着传入的MS2即没有PA也没有FA
  havePA <- length(MS2$PA)
  haveFA <- length(MS2$FA)
  #计算MS2是没有PA还是没有FA，还是两者都有
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if ((havePA == 1) & (haveFA == 5)) {
    spliceFA <- purrr::pmap(.l = list(as.list(MS1$`Chain Length:Δ`) , as.list(MS1$Oxform)) , .f = findmyCL::splice4FA , FA = MS2$FA)
  #进行4个FA的拼接
    spliceFA <- findmyCL::deleteNULL(spliceFA)
    #去除NULL
    if (length(spliceFA) == 0) {
      return(NULL)
    }
    #如果没有拼接成功，则返回空值
    #如果拼接成功，进行下面步骤
    names(spliceFA) <- 1:length(spliceFA)
    #重命名
    MS2 <- findmyCL::fappend(list1 = MS2 , list2 = spliceFA)
    #追加到MS2中
    return(MS2)
  }
  #假如没有PA，则只拼接FA,即只要4个FA能组成一个心磷脂就好
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if ((haveFA == 1) & (havePA == 5)) {
    splicePA_result <- purrr::pmap(.l = list(as.list(MS1$`Chain Length:Δ`) , as.list(MS1$Oxform)) , .f = findmyCL::splice2PA , PA = MS2$PA)
    #进行PA的拼接
    splicePA_result <- findmyCL::deleteNULL(splicePA_result)
    #去除NULL
    if (length(splicePA_result) == 0) {
      return(NULL)
    }
    #如果这个二级的PA无法拼接成心磷脂，则将这个二级去除（即返回NULL）
    #如果拼接成功，进行下面步骤
    names(splicePA_result) <- 1:length(splicePA_result)
    #重命名
    splicePA <- splicePA_result
    #赋值到splicePA方便命名
    MS2 <- findmyCL::fappend(list1 = MS2 , list2 = splicePA)
    #追加到MS2中
    return(MS2)
  }
  #假如没有FA，则只拼接PA
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if ((haveFA == 5) & (havePA == 5)) {
    splicePA_result <- purrr::pmap(.l = list(as.list(MS1$`Chain Length:Δ`) , as.list(MS1$Oxform)) , .f = findmyCL::splice2PA , PA = MS2$PA)
    #进行PA的拼接
    splicePA_result <- findmyCL::deleteNULL(splicePA_result)
    #去除NULL
    if (length(splicePA_result) == 0) {
      return(NULL)
    }
    #如果这个二级的PA都无法拼接成心磷脂，则将这个二级去除（即返回NULL）
    names(splicePA_result) <- 1:length(splicePA_result)
    #重命名拼接PA的结果

    #future::plan(multisession , workers = 4)
    spliceFA <- purrr::map(.x = splicePA_result , .f = findmyCL::splice2FA , MS2 = MS2)
    #将每个splicePA_result的子集计算所有FA
    spliceFA <- findmyCL::deleteNULL(list = spliceFA)
    #去除空值
    if (length(spliceFA) == 0) {
      return(NULL)
    }
    #如果FA的拼接没有结果（两条PA都没有拼接出FA或者有一条没有拼接出FA），则返回空值

    splicePA <- splicePA_result
    #赋值到splicePA中免得命名

    MS2 <- findmyCL::fappend(list1 = MS2 , list2 = splicePA)
    MS2 <- findmyCL::fappend(list1 = MS2 , list2 = spliceFA)
    #追加到MS2中

    return(MS2)
    #返回结果
  }
  #假如既有PA，也有FA，则先拼接PA，再拼接FA
}
#三种可能，只有PA，只有FA，既有FA又有PA

#' @title Splice 2PA into a cardiolipin
#' @param PA dataframe(1)
#' @param chain_double character(1)
#' @param oxygen numeric(1)
#' @return list(1),including splicePA-result,its corresponding Chain Length:Δ && oxygen
#' @export
splice2PA <- function(PA , chain_double , oxygen){
  PA_row_length <- length(PA[,1])
  #求有多少种PA
  sample_number <- 2
  #任意两个相加
  factor <- gl(PA_row_length , sample_number) %>%
    as.vector() %>%
    as.numeric()
  #生成行号的可能的组合情况
  combination <- combn(factor , m = sample_number) %>%
    t() %>%
    unique() %>%
    t()
  #组合，并除去重复的情况
  add_chain_result <- add_2Chain(vector1 = PA[combination[1,],4] , vector2 = PA[combination[2,],4])
  #将两个向量的“a:b”“c:d”相加，成为“a+c：c+d”
  add_oxygen_result <- PA[combination[1,],5] + PA[combination[2,],5]
  #将对于的氧原子个数相加
  matchresult <- intersect(x = which(add_chain_result == chain_double) , y = which(add_oxygen_result == oxygen))
  #即满足chain相加能与一级的相同，即满足PA的氧原子个数相加能与一级的氧原子个数相同的行号
  if (length(matchresult) == 0) {
    return(NULL)
  }
  #如果PA没有拼接成功，则返回NULL
  if (length(matchresult) > 0) {
    store <- list()
    #存储最后结果
    matchresult_message <- combination[ , matchresult]
    #将配对成功的行号保存
    matchresult_message <- as.matrix(matchresult_message)
    #转变为矩阵方便后续操作
    rbind_result <- purrr::pmap(.l = list(as.list(matchresult_message[1,]) , as.list(matchresult_message[2,])) , .f = findmyCL::catchDataframe2Col , dataframe = PA) %>%
      rbind.fill()
    #将配对成功的对于行号的行绑定在一起,组合成大数据框
    colnames(rbind_result) <- rep(colnames(PA),2)
    #重命名结果数据框

    store <- findmyCL::fappend(list1 = store,list2 = rbind_result)
    #追加到结果中
    store <- append(x = store , values = chain_double)
    store <- append(x = store , values = oxygen)
    #追加chian:double和氧原子到结果中

    names(store) <- c("PA" , "Chain Length:Δ" , "oxygen")
    #重命名结果
    return(store)
  }
  #如果拼接成功，返回结果
}
#判断PA的结果中哪些可以拼接成心磷脂,如果PA不能拼接不成心磷脂，返回结果0

#' @title Bind 2 row by column in a dataframe
#' @param row1 numeric(1)
#' @param row2 numeric(1)
#' @param dataframe dataframe(1)
#' @return dataframe(1)
#' @export
catchDataframe2Col <- function(row1 , row2 , dataframe){
  dataframe <- cbind(dataframe[row1,] , dataframe[row2,])
  colnames(dataframe) <- 1:length(colnames(dataframe))
  #重命名数据框方便后续合并多个数据框

  return(dataframe)
}
#绑定一个数据框的某两行

#' @title Bind 4 row by column in a dataframe
#' @param row1 numeric(1)
#' @param row2 numeric(2)
#' @param row3 numeric(3)
#' @param row4 numeric(4)
#' @param dataframe dataframe(1)
#' @return dataframe(1)
#' @export
catchDataframe4Col <- function(row1 , row2 , row3 , row4 ,dataframe){
  dataframe <- cbind(dataframe[row1,] , dataframe[row2,] , dataframe[row3,] , dataframe[row4,])
  #合并4行为一行
  return(dataframe)
}
#绑定一个数据框的某四行

#' @title Splice 2 FA from the result in PA-splicing
#' @param splicePA list(1)
#' @param MS2 list(1)
#' @export
splice2FA <- function(splicePA , MS2){
  spliceFA_result <- splice2FA_main(splicePA_dataframe = splicePA[[1]] , MS2 = MS2)
  #计算匹配到的FA，输入的是数据框和二级
  if (is.null(spliceFA_result)) {
    return(NULL)
  }
  #如果一种PA没有匹配到FA，则返回NULL
  spliceFA_result <- findmyCL::deleteNULL(spliceFA_result)
  #去除NULL
  if (length(spliceFA_result) == 0) {
    return(NULL)
  }
  #如果FA没有配对成功，则返回空值

  if (length(spliceFA_result) > 0) {
    return(spliceFA_result)
    #返回结果
  }
}

#' @title Splice 2 FA from the result in PA-splicing
#' @param splicePA_dataframe dataframe(1)
#' @param MS2 list(1)
#' @export
splice2FA_main <- function(splicePA_dataframe , MS2){
  data <- t(splicePA_dataframe)
  vector <- lapply(seq_len(ncol(data)), function(i) data[,i])
  #将数据框转换为列表便于多线程
  spliceFA <- purrr::map(.x = vector , .f = findmyCL::splice2FA_vector , MS2 = MS2)
  #对数据框中的每一行进行FA的匹配
  spliceFA <- findmyCL::deleteNULL(spliceFA)
  #去除NULL
  if (length(spliceFA) == 0) {
    return(NULL)
  }
  #如果没有FA配对成功 or 两条PA只有一条PA能拼出FA，返回空值
  if (length(spliceFA) > 0) {
    chain_double <- findmyCL::add_2Chain(splicePA_dataframe[1,4] , splicePA_dataframe[1,9])
    oxygen <- splicePA_dataframe[1,5] + splicePA_dataframe[1,10]
    #计算chain_double和氧
    names(spliceFA) <- 1:length(spliceFA)
    #重命名(不加后面会有bug)
    spliceFA <- append(spliceFA , values = c(chain_double , oxygen))
    #追加chain_double和oxygen
    names(spliceFA)[(length(spliceFA) - 1):(length(spliceFA))] <- c("Chain Length:Δ" , "oxygen")
    #重命名
    return(spliceFA)
  }
  #如果配对FA成功，返回配对结果
}
#配对单个PA配对结果的数据框是否有FA,输入的是PA配对之后的数据框，成对出现， 如果没有FA配对成功 or 两条PA只有一条PA能拼出FA，返回空值

#' @title Splice 2 FA from the result in PA-splicing
#' @param vector vector(1)
#' @param MS2 list(1)
#' @export
splice2FA_vector <- function(vector , MS2){
  result1 <- findmyCL::splice2PA(PA = MS2$FA , chain_double = vector[4]  , oxygen = vector[5])
  result2 <-  findmyCL::splice2PA(PA = MS2$FA , chain_double = vector[9]  , oxygen = vector[10])
  if ((length(result1) == 0)||(length(result2) == 0)) {
    return(NULL)
  }
  #如果任意一条PA没有拼接成FA，则将这两条PA去除
  if ((length(result1) > 0) && (length(result2) > 0)) {
    result_colname <- colnames(result1[["PA"]])
    colnames(result1[["PA"]]) <- 1:length(colnames(result1[["PA"]]))
    colnames(result2[["PA"]]) <- 11 : (10 + (length(colnames(result2[["PA"]]))))
    if ((vector[4] != vector[9]) || (vector[5] != vector[10]) ) {
      combine_FA <- tidyr::crossing(result1[["PA"]] , result2[["PA"]])
      #组合FA
    }
    #如果两个PA的chain：▲不一样或者氧原子个数不一样（就是两个PA是不一样的PA），则进行FA的任意组合
    if ((vector[4] == vector[9]) && (vector[5] == vector[10])) {
      help_dataframe1 <- 1:length(result1[["PA"]][,1]) %>%
        as.data.frame()
      colnames(help_dataframe1) <- 1

      help_dataframe2 <- 1:length(result1[["PA"]][,1]) %>%
        as.data.frame()
      colnames(help_dataframe2) <- 2

      help_dataframe <- tidyr::crossing(help_dataframe1,help_dataframe2) %>%
        apply(1 , sort , decreasing = T) %>%
        t() %>%
        as.data.frame() %>%
        dplyr::distinct()
      #创建帮助组合PA相同的数据框（里面是行号的排列组合）

      dataframe1 <- result1[["PA"]][help_dataframe[,1] , ]
      dataframe2 <- result2[["PA"]][help_dataframe[,2] , ]
      combine_FA <- cbind(dataframe1 , dataframe2)
      #将两个结果中对应的行行组合在一起，形成一个大数据框
    }
    #如果两个PA的chain：▲和氧原子个数都一，则进行FA的任意组合（重复的去除）
    colnames(combine_FA) <- rep(result_colname , 2)
    #重命名结果
    as.data.frame(combine_FA) %>%
     return()
  }
  #如果都能拼接成FA，则组合所有FA
}
#输入一个数据框中的向量(就是配对成功的两条PA放在同一行中了，需要这一行PA)，返回可能的FA，已经进行了两对PA的FA的任意组合

#' @title Splice 4 FA into a cardiolipin
#' @param FA dataframe(1)
#' @param chain_double character(1)
#' @param oxygen numeric(1)
#' @return list(1)
#' @export
splice4FA <- function(FA , chain_double , oxygen){
  FA_num <- length(FA[,1])
  #求一共有多少个FA
  FA_combination <- expand.grid(rep(list(1 : FA_num), 4)) %>%
    apply(1 , sort , decreasing = T) %>%
    t() %>%
    as.data.frame() %>%
    dplyr::distinct()
  #计算所有的组合可能

  chain_add_result <- findmyCL::add_2Chain(vector1 = FA[FA_combination[ , 1] , 4] , FA[FA_combination[ , 2] , 4]) %>%
    findmyCL::add_2Chain(FA[FA_combination[ , 3] , 4]) %>%
    findmyCL::add_2Chain(FA[FA_combination[ , 4] , 4])
  #将所有组合可能的chain:double相加，4行相加
  oxygen_add_result <- FA[FA_combination[ , 1] , 5] + FA[FA_combination[ , 2] , 5] + FA[FA_combination[ , 3] , 5] + FA[FA_combination[ , 4] , 5]
  #将所有组合可能的氧原子个数相加，4行相加
  matchresult <- which(chain_add_result == chain_double) && which(oxygen_add_result == 14)
  #即与chain：double相等又与氧原子个数相等的行号
  if (is.na(matchresult)) {
    return(NULL)
  }
  #如果配对没有成功，则返回空值

  FA_combination_result <- FA_combination[matchresult , ]
  #将配对成功结果的行号保存
  catch_result <- purrr::pmap(.l = list(row1 = list(FA_combination_result[,1]) , row2 = list(FA_combination_result[,2]) , row3 = list(FA_combination_result[,3]) , row4 = list(FA_combination_result[,4])) , .f = findmyCL::catchDataframe4Col , dataframe = FA)
  #将所有符合的行提取
  catch_result <- append(x = catch_result , values = chain_double)
  catch_result <- append(x = catch_result , values = oxygen)
  names(catch_result) <- c("FA","Chain Length:Δ" , "oxygen")

  return(catch_result)
  #
}
#如果没有PA但是有PA，则只拼接4条FA
