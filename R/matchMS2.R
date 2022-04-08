#' @title Matching the MS2 database
#' @description This function use result in MS2-checkresult to search PA&FA and splice Cardiolipin.The MS2 database consists of PA&FA database.Before searching PA&FA,we use  [findmyCL::check_MS1_MS2_HeadFAPA()] to check if a MS2 have a head group(mz:152.9963485799).If this MS2 has not head group,we would delete it.If this MS2 have head group,then we would find PA&FA in this MS2.After searching PA&FA,we would splice PA&FA into Cardiolipin [findmyCL::splicePAFA()] .If this MS2 cannot splice into Cardiolipin,we would delete this MS2.Moreover,we would delete a MS1 and its MS2 if we connot find Cardiolipin in all its MS2.
#' @param object a findCL object after [findmyCL::checkMS2()]
#' @param ppm numeric(1) defining the maximal tolerated m/z deviation in consecutive scans in parts per million (ppm) for the initial ROI definition when match our database.Default is 30.
#' @param workers numeric(1) Number of multithreaded cores you use in this function.Default is 2.
#' @import purrr
#' @import furrr
#' @import future
#' @import dplyr
#' @import progressr
#' @return
#' @export
matchMS2 <- function(object , ppm = 30 , workers = 2){
  ms2checkresult_name <- as.list(names(object@ms2CheckResult))
  ms2checkresult_name <- head(ms2checkresult_name , -1)
  #去除ppm属性
  # future::plan(multisession , workers = 4)
  # ms2MatchResult <- furrr::future_map(.x = ms2checkresult_name , .f = findmyCL::matchMS2_main , object = object , ppm = ppm)
  ms2MatchResult <- purrr::map(.x = ms2checkresult_name , .f = findmyCL::matchMS2_main , object = object , ppm = ppm , workers = 2)
  #计算所有ms2CheckResult的结果
  names(ms2MatchResult) <- ms2checkresult_name
  #重命名
  ms2MatchResult <- findmyCL::deleteNULL(ms2MatchResult)
  #去除NULL
  ms2MatchResult <- append(ms2MatchResult,ppm)
  names(ms2MatchResult)[length(ms2MatchResult)] <- "ppm"
  #加入ppm属性
  object@ms2MatchResult <- ms2MatchResult
  #加入到结果中
  rtfilter_result <- findmyCL::rtfilter(object = object)
  #把对应有多个心磷脂的一级和只对应一个心磷脂的一级分开，并且把多余的splicePA的结果去除，1对1的作为种子，1对多的要被判断
  object@ms2MatchResult <- rtfilter_result
  #将最终结果加入object中
  return(object)
}

#' @title Search every MS2-checking for their PA&FA and splice them.
#' @description search every MS2-checking for their PA&FA and splice them,every element use [check_MS1_MS2_HeadFAPA()] to search and splice.
#' @param object findmyCLclass(1)
#' @param ppm numeric(1)
#' @param ms2checkresult_name character(1)
#' @param workers numeric(1) Number of cores used by multithreading
#' @export
matchMS2_main <- function(object , ppm = 30 ,ms2checkresult_name = "CL" ,workers = 2){
  dealing_list <- object@ms2CheckResult[[ms2checkresult_name]]
  #将待处理的列表输出
  paste0("Searing and splicing from: ",ms2checkresult_name) %>%
    message()
  #提示信息
  future::plan(multisession , workers = workers)
  #设置多线程使用核数

  with_progress({
    p <- progressr::progressor(steps = length(dealing_list))
    #设置进度条
    spliceResult <- furrr::future_map(.x = dealing_list , .f = findmyCL::check_MS1_MS2_HeadFAPA , ppm = ppm , ms2checkresult_name = ms2checkresult_name , p = p)
    #迭代输出结果,检查了头基并且拼接了心磷脂
  })

  cat('\n')
  #换行
  message("Done!")
  #提示信息(完成了这一组心磷脂的拼接)
  spliceResult <- findmyCL::deleteNULL(list = spliceResult)
  #去除空值
  if (length(spliceResult) == 0) {
    return(NULL)
  }
  #如果都没有头基，则直接将这个结果删除
  names(spliceResult) <- 1:length(spliceResult)
  #重命名结果
  return(spliceResult)
}
#配对ms2checkresult中比对的每个结果,比如“CL”结果,然后对结果进行拼接成心磷脂。spliceResult_esr

#' @title Checking if a big MS1 have its corresponding MS2
#' @description checking if a aggregation which have MS1 and MS2 have PA or FA and splice PA and FA in cardiolipin.
#' @param MS1 list(1)
#' @param ppm numeric(1)
#' @param p progress bar
#' @export
check_MS1_MS2_HeadFAPA <- function(bigMS1, ppm = 30 , ms2checkresult_name , p){
  p()
  #显示进度条
  MS2_dealing <- bigMS1[-1]
  #去除一级的信息，剩下的是二级的信息
  result <- purrr::map(.x = MS2_dealing,.f = findmyCL::check_HeadFAPA,ppm = ppm,MS1 = bigMS1$MS1 , ms2checkresult_name = ms2checkresult_name)
  #批量检查所有一个MS1中的所有MS2是否有头基和PA、FA，如果没有头基，则返回NULL，如果有头基但是没有PA和FA，也返回NULL

  result <- findmyCL::deleteNULL(list = result)
  #去除空值(NULL)
  if (length(result) == 0) {
    return(NULL)
  }
  #如果这些二级都没有头基，则将这个一级和这些二级一起删除，返回空值
  if (length(result) > 0) {
    store <- list()
    #存储结果
    store <- findmyCL::fappend(store,bigMS1[[1]])
    #追加一级的信息
    store <- append(x = store,values = result)
    #存储结果
    names(store)[1] <- "MS1"
    return(store)
  }
}
#检查一个包含一个一级和多个二级的集合的二级是否有头基，如果二级都没有头基，则将这个集合去除

#' @title Checking PA and FA for a MS2
#' @description Checking if a MS2 have its corresponding PA or FA and splice them,use [findPAFA()] to find PA and FA.
#' @param MS1 list(1)
#' @param MS2 list(1)
#' @param ppm numeric(1)
#' @seealso [splicePAFA()]
#' @export
check_HeadFAPA <- function(MS1 , MS2 , ppm = 30 , ms2checkresult_name){
  min_match_1 <- 152.9963485799 - 152.9963485799*ppm/1000000
  max_match_1 <- 152.9963485799 + 152.9963485799*ppm/1000000
  #求头基出现的最大mz、最小mz,152.9963485799为精确头基mz,一共有1种
  match_result_1 <- length(which((MS2@mz > min_match_1) & (MS2@mz < max_match_1)))
  match_result <- match_result_1
  #计算这个二级中是否有心磷脂的头基

  if (match_result == 0) {
    return(NULL)
  }
  #如果没有头基，返回空值

  if (match_result > 0) {
    MS2 <- findmyCL::findPAFA(MS2 = MS2 , ppm = ppm)
    #如果有头基，则寻找PA和FA，这里算出的MS2不会同时没有PA和FA,如果同时没有PA和FA的话，那么MS2就会是NULL
    if (ms2checkresult_name == "CL") {
      MS2 <- findmyCL::splicePAFA(MS2 = MS2 , MS1 = MS1)
      #拼接FA、PA，如果没有拼接成功，那么MS2就会是NULL
      return(MS2)
    }
    #如果是CL（即有4条FA，则按照splicePAFA的方法来拼接成心磷脂）
    if (ms2checkresult_name == "MLCL") {
      MS2 <- findmyCL::splicePAFA_MLCL(MS2 = MS2 , MS1 = MS1)
      #拼接FA、PA，如果没有拼接成功，那么MS2就会是NULL
      return(MS2)
    }
    #如果是MLCL（即有3条FA，则按照splicePAFA_MLCL的方法来拼接成心磷脂）
    if (ms2checkresult_name == "DLCL") {
      MS2 <- findmyCL::splicePAFA_DLCL(MS2 = MS2 , MS1 = MS1)
      return(MS2)
    }
    #如果是DLCL（即有2条FA，则按照splicePAFA_DLCL的方法来拼接成心磷脂）
  }
  #如果有头基，则配对PA和FA，如果没有PA和FA，则返回结果MS2为NULL
}
#检查一个二级是否有头基,如果没有头基，则将其删除，如果有头基，没有PA和FA，也返回NULL,并且进行了心磷脂的拼接，拼接可看splicePAFA

#' @title Finding PA and FA in a MS2
#' @param MS2 list(1)
#' @param ppm numeric(1)
#' @export
findPAFA <- function(MS2 , ppm = 30){
  MS2_dataframe <- cbind(MS2@mz,MS2@intensity) %>%
    as.data.frame()
  #创建用于比对的数据
  colnames(MS2_dataframe) <- c("mz","intensity")
  #改列名
  MS2_dataframe <- t(MS2_dataframe)
  vector <- lapply(seq_len(ncol(MS2_dataframe)), function(i) MS2_dataframe[,i])
  #将数据框转换为列表便于多线程
  PA <- purrr::map(.x = vector , .f = findmyCL::matchDatabaseMS2 , vector_column = 1 , ppm = ppm , database = findmyCL:::MS2_database$PA , database_column = 1)
  #迭代计算配对的PA
  FA <- purrr::map(.x = vector , .f = findmyCL::matchDatabaseMS2 , vector_column = 1 , ppm = ppm , database = findmyCL:::MS2_database$FA , database_column = 1)
  #迭代计算配对的FA,使用函数matchDatabaseMS2是为了不引入进度条

  PA <- findmyCL::deleteNULL(list = PA)
  FA <- findmyCL::deleteNULL(list = FA)
  #去除NULL
  if ((length(PA) == 0) & (length(FA) == 0)) {
    return(NULL)
  }
  #如果PA和FA都没有配对上的值，则返回空值
  if (length(PA) != 0) {
    PA <- plyr::rbind.fill(PA)
    colnames_length <- length(rownames(MS2_dataframe))
    #计算原本数据列名长度
    result_name_length <- length(colnames(PA))
    #计算结果的列名长度
    colnames(PA)[(colnames_length + 1):result_name_length] <- c("aditive form" , "Chain Length:Δ" , "Oxform")
    colnames(PA)[1:colnames_length] <- rownames(MS2_dataframe)
    #改结果的列名
    PA_chain <- PA$`Chain Length:Δ` %>%
      strsplit(":") %>%
      as.data.frame() %>%
      t() %>%
      as.data.frame()
    #提取PA的碳数与不饱和度
    colnames(PA_chain) <- c("Chain Length" , "Δ")
    #重命名提取的PA的碳数与不饱和度
    oxygen_limit <- ((as.numeric(PA_chain$`Chain Length`) - as.numeric(PA_chain$Δ)*2 - 2) / 2) %>%
      ceiling()
    #计算允许的最大氧的个数
    PA_row <- which(PA$Oxform <= oxygen_limit)
    #能够符合规律的FA的行数
    if (length(PA_row) != 0) {
      PA <- PA[PA_row ,]
      #重新赋值到FA
    }
    #如果存在不符合规律的PA，则将其去除
    rownames(PA) <- 1:length(PA[,1])
    #重命名PA的行名
  }
  #用限制条件限制PA的氧个数
  if (length(PA) == 0) {
    PA <- 0
  }
  if (length(FA) != 0) {
    FA <- plyr::rbind.fill(FA)
    colnames_length <- length(rownames(MS2_dataframe))
    #计算原本数据列名长度
    result_name_length <- length(colnames(FA))
    #计算结果的列名长度
    colnames(FA)[(colnames_length + 1):result_name_length] <- c("aditive form" , "Chain Length:Δ" , "Oxform")
    colnames(FA)[1:colnames_length] <- rownames(MS2_dataframe)
    #改结果的列名
    FA_chain <- FA$`Chain Length:Δ` %>%
      strsplit(":") %>%
      as.data.frame() %>%
      t() %>%
      as.data.frame()
    #提取FA的碳数与不饱和度
    colnames(FA_chain) <- c("Chain Length" , "Δ")
    #重命名提取的FA的碳数与不饱和度
    oxygen_limit <- ((as.numeric(FA_chain$`Chain Length`) - as.numeric(FA_chain$Δ)*2 - 1) / 2) %>%
      ceiling()
    #计算允许的最大氧的个数
    FA_row <- which(FA$Oxform <= oxygen_limit)
    #能够符合规律的FA的行数
    if (length(FA_row) != 0) {
      FA <- FA[FA_row ,]
      #重新赋值到FA
    }
    #如果存在不符合规律的FA，则将其去除

    rownames(FA) <- 1:length(FA[,1])
    #重命名FA的行名
  }
  #用限制条件限制FA的氧个数
  if (length(FA) == 0) {
    FA <- 0
  }
  #处理PA和FA
  result <- findmyCL::fappend(list1 = MS2,list2 = PA)
  result <- findmyCL::fappend(list1 = result,list2 = FA)
  #将结果加入到result中
  names(result)[1] <- "MS2"
  #重命名结果
  return(result)
}
#配对一个MS2中的PA和FA

