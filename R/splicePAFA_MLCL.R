#' @title Splice PA and FA from a MS2 and MS1 in MLCL
#' @description A function to search PA and FA from MS2.If this MS1 can not search its 2 PA,delete this MS1(this MS1 may not have just one corresponding Cardiolipin).If this MS2 can search its PA but can not searth corresponding 2 FA or it can just match one FA,delete this MS2.
#' @param MS2 list(1),including MS2 information,PA&FA search result
#' @param MS1 dataframe(1),including mz,rt,rtmin,rtmax ...
#' @return list(1),a new MS2,adding PA&FA splicing result
#' @export
splicePAFA_MLCL <- function(MS2 , MS1){
  if (is.null(MS2)) {
    return(NULL)
  }
  #如果传入的MS2是空值，则直接返回空值，因为这代表着传入的MS2即没有PA也没有FA
  havePA <- length(MS2$PA)
  haveFA <- length(MS2$FA)
  #计算MS2是没有PA还是没有FA，还是两者都有
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if ((havePA == 1) & (haveFA == 5)){

  }
  #如果没有PA，则只拼接FA,即只要3个FA能组成一个MLCL就好
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if ((haveFA == 1) & (havePA == 5)){
    return(NULL)
  }
  #假如没有FA，则只拼接PA，但是拼接成MLCL需要一个PA+一个FA，若没有FA，则无法拼接成MLCL
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if ((haveFA == 5) & (havePA == 5)){
    splicePA_result <- purrr::pmap(.l = list(as.list(MS1$`Chain Length:Δ`) , as.list(MS1$Oxform)) , .f = findmyCL::splicePA_MLCL , PA = MS2$PA , FA = MS2$FA)
    splicePA_result <- findmyCL::deleteNULL(splicePA_result)
    #去除NULL
    if (length(splicePA_result) == 0) {
      return(NULL)
    }
    #如果这个二级的PA都无法拼接成心磷脂，则将这个二级去除（即返回NULL）
    names(splicePA_result) <- 1:length(splicePA_result)
    #重命名拼接PA的结果
    spliceFA_result <- purrr::map(.x = splicePA_result , .f = findmyCL::splice2FA_MLCL , FA = MS2$FA) %>%
      findmyCL::deleteNULL()
    #用PA的拼接结果拼接2个FA，并且找另外一个FA能一起拼接成MLCL
    if (length(spliceFA_result) == 0) {
      return(NULL)
    }
    #如果不能拼接成MLCL，则返回NULL
    return(spliceFA_result)
    #返回最终结果
  }
  #假如既有PA，也有FA，则先拼接成一个PA，PA再由FA拼接,再找1个FA拼接成MLCL
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
}

#' @title Search which PA can splice into a MLCL with a FA.
#' @description Search which PA can splice into a MLCL with a FA.
#' @param PA dataframe(1)
#' @param chain_double character(1)
#' @param oxygen numeric(1)
#' @return list(1),including splicePA-result,its corresponding Chain Length:Δ && oxygen
#' @export
splicePA_MLCL <- function(PA , chain_double , oxygen , FA){
  PA_length <- length(PA[,1])
  #求有多少个PA
  store <- c()
  #存储不需要的PA
  for (i in 1:PA_length) {
    chain_add_result <- findmyCL::add_2Chain(vector1 = PA[i,]$`Chain Length:Δ` , vector2 = FA$`Chain Length:Δ`)
    #PA中的Chain Length:Δ与FA中的Chain Length:Δ相加
    oxygen_add_result <- PA[i,]$Oxform + FA$Oxform
    #PA中的Oxform和FA中的Oxform相加
    chain_matchresult <- which(chain_add_result == chain_double)
    oxygen_matchresult <- which(oxygen_add_result == oxygen)
    #求相加的结果与一级中的Chain Length:Δ、oxygen相符合的结果
    if ((length(chain_matchresult) == 0) || (length(oxygen_matchresult) == 0)) {
      store <- c(store , i)
    }
    #如果相加结果中的Chain Length:Δ或者oxygen无法和一级的配对上，则存储再store中
  }

  if (length(store) == PA_length) {
    return(NULL)
  }
  #如果PA和FA无法拼接成MLCL，则返回NULL

  PA <- PA[-store,]
  rownames(PA) <- 1:length(PA[,1])
  #去除无法拼接成MLCL的PA，重命名行名

  result_list <- list(PA , chain_double , oxygen)
  names(result_list) <- c("PA" , "Chain Length:Δ" , "oxygen")
  #重命名结果
  return(result_list)
  #返回结果
}

#' @title Splice 1PA & 1FA into MLCL.
#' @description First,splice 2FA into 1PA.Second,search the third FA which can be splice into a MLCL with this 2FA(can be splice into 1PA).
#' @param splicePA list(1),including PA,Chain Length:Δ,oxygen.
#' @param FA dataframe(1)
#' @return list(1)
#' @export
splice2FA_MLCL <- function(splicePA , FA){
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
  spliceFA_result <- plyr::rbind.fill(spliceFA_result)
  #从列表转换为数据框，并且将所有结果行整合

  list <- findmyCL::turnDataframeList(dataframe = spliceFA_result)
  #将数据框转换为列表便于多线程
  spliceFA_result <- purrr::map(.x = list , .f = findmyCL::splice1FA_1PA_into_MLCL , FA = FA ,chain_double = splicePA$`Chain Length:Δ` , oxygen = splicePA$oxygen) %>%
    findmyCL::deleteNULL()
  #计算2个FA和1个FA能否组成一个MLCL
  if (length(spliceFA_result) == 0) {
    return(NULL)
  }
  #如果无法拼接成MLCL，则返回NULL
  spliceFA_result <- plyr::rbind.fill(spliceFA_result)
  #从列表转换为数据框，并且将所有结果行整合
  store <- list(spliceFA_result , splicePA[["Chain Length:Δ"]] , splicePA[["oxygen"]])
  #存储FA拼接的结果
  names(store) <- c("FA" , "Chain Length:Δ" , "oxygen")
  #重命名结果
  return(store)
}

#' @title Splice 2FA into a PA in MLCL
#' @description Use [findmyCL::splice2PA()] to splice 2FA into a PA,return a dataframe.
#' @param splicePA_dataframe dataframe(1),have 1 row.
#' @param FA dataframe(1)
#' @return dataframe(1)
#' @export
splice2FA_MLCL_main <- function(splicePA_dataframe , FA){
  if (class(splicePA_dataframe) == "character") {
    splicePA_dataframe <- as.data.frame(splicePA_dataframe) %>%
      t() %>%
      as.data.frame()
  }
  #如果输入的是一个向量，则需要先转换为数据框
  result <- findmyCL::splice2PA(PA = FA , chain_double = splicePA_dataframe$`Chain Length:Δ` , oxygen = splicePA_dataframe$Oxform)
  #利用splice2PA函数搜索这个PA能由FA中的哪两个FA组成
  spliceFA <- result[["PA"]]
  #去除拼接的结果，结果位于PA的槽中
  return(spliceFA)
  #返回结果
}

#' @title Splice 1FA and 1PA into a MLCL
#' @description Give a PA which is spliced by 2FA.Use this PA to search which FA in FA(dataframe) can be splice into a MLCL together.
#' @param spliceFA_result_dataframe dataframe(1),have 1 row
#' @param FA dataframe(1)
#' @param chain_double character(1),Chain Length:Δ from MS1
#' @param oxygen numeric(1),Oxform from MS1
#' @return dataframe(1)
#' @export
splice1FA_1PA_into_MLCL <- function(spliceFA_result_dataframe , FA , chain_double , oxygen){
  spliceFA_result_dataframe <- as.matrix(spliceFA_result_dataframe) %>%
    t() %>%
    as.data.frame()
  #转换为数据框，因为传入的可能是向量
  chain_add_result <- findmyCL::add_2Chain(vector1 = spliceFA_result_dataframe[1,4] , vector2 = spliceFA_result_dataframe[1,9]) %>%
    findmyCL::add_2Chain(vector2 = FA$`Chain Length:Δ`)
  #将2个FA的Chain Length:Δ和单个FA的Chain Length:Δ相加
  Oxform_add_result <- as.numeric(spliceFA_result_dataframe[1,5]) + as.numeric(spliceFA_result_dataframe[1,10]) + FA$Oxform
  #将2个FA的Oxform之和与单个FA的Oxform相加
  match_number <- which((chain_add_result == chain_double) && (Oxform_add_result == oxygen))
  #求能与一级的Chain Length:Δ和Oxform相符的结果
  if (length(match_number) == 0) {
    return(NULL)
  }
  #如果3个FAChain Length:Δ相加和Oxform相加的结合没有与一级心磷脂的Chain Length:Δ和Oxform对得上，则返回空值
  options(warn = -1)
  #不提示警告
  result <- cbind(spliceFA_result_dataframe , FA[match_number,])
  options(warn = 1)
  #提示警告
  return(result)
  #返回结果
}


#' @title Splice 3FA into a MLCL
#' @description Splice 3FA into a MLCL which have no PA.
#' @param FA dataframe(1)
#' @param chain_double character(1)，like "14:0".
#' @param oxygen numeric(1)
#' @return list(1)
#' @export
splice3FA_MLCL <- function(FA , chain_double , oxygen){
  FA_num <- length(FA[,1])
  #求一共有多少个FA
  FA_combination <- expand.grid(rep(list(1 : FA_num), 3)) %>%
    apply(1 , sort , decreasing = T) %>%
    t() %>%
    as.data.frame() %>%
    dplyr::distinct()
  #计算所有的组合可能(行号的组合)
  chain_add_result <- findmyCL::add_2Chain(vector1 = FA[FA_combination[ , 1] , 4] , FA[FA_combination[ , 2] , 4]) %>%
    findmyCL::add_2Chain(FA[FA_combination[ , 3] , 4])
  #将所有组合可能的chain:double相加，3行相加(FA中的chain:double相加)
  oxygen_add_result <- FA[FA_combination[ , 1] , 5] + FA[FA_combination[ , 2] , 5] + FA[FA_combination[ , 3] , 5]
  #将所有组合可能的氧原子个数相加，3行相加(FA中的oxygen相加)
  matchresult <- intersect(which(chain_add_result == chain_double) , which(oxygen_add_result == oxygen))
  #取即与chain：double相等又与氧原子个数相等的行号
  if (length(matchresult) == 0) {
    return(NULL)
  }
  #如果配对没有成功，则返回空值
  FA_combination_result <- FA_combination[matchresult , ]
  #将配对成功结果的行号保存
  catch_result <- purrr::pmap(.l = list(row1 = list(FA_combination_result[,1]) , row2 = list(FA_combination_result[,2]) , row3 = list(FA_combination_result[,3])) , .f = findmyCL::catchDataframe3Col , dataframe = FA)
  #将所有符合的行提取
  catch_result <- append(x = catch_result , values = chain_double)
  catch_result <- append(x = catch_result , values = oxygen)
  names(catch_result) <- c("FA","Chain Length:Δ" , "oxygen")
  #追加属性，重命名结果
  return(catch_result)
  #返回结果
}

#' @title Bind 3 row by column in a dataframe
#' @param row1 numeric(1)
#' @param row2 numeric(2)
#' @param row3 numeric(3)
#' @param row4 numeric(4)
#' @param dataframe dataframe(1)
#' @return dataframe(1)
#' @export
catchDataframe3Col <- function(row1 , row2 , row3 ,dataframe){
  dataframe <- cbind(dataframe[row1,] , dataframe[row2,] , dataframe[row3,] )
  #合并4行为一行
  return(dataframe)
}
#绑定一个数据框的某四行
