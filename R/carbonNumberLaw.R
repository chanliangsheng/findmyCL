#' @title Use carbon number law to filter some wrong Cardiolipin.
#' @param dataframe dataframe(1) including chain,Δ,rt like (56:0,5).Isomer is allowed.
#' @return dataframe(1) including chain,Δ,rt like (56:0,5) after filtering.
#' @export
carbonNumberLaw <- function(dataframe) {
  for (i in 1:length(dataframe[, 1])) {
    score <- findmyCL::score(dataframe = dataframe)
    #对心磷脂数据进行打分
    if (length(which(score$score != 1)) == 0) {
      return(dataframe)
    }
    #如果所有心磷脂的打分都为1，则直接返回这组心磷脂
    if (length(which(score$score != 1)) > 0) {
      delete_row <- which(score$score == min(score$score), arr.ind = TRUE)

      if (length(delete_row) > 1) {
        delete_row <- sample(x = delete_row , size = 1)
      }
      #查找分数最小的心磷脂的行(可能有多个分数一样小的心磷脂),如果有多个分数一样小的，则随机选取一个进行去除
      dataframe <- dataframe[-delete_row , ]
      #删除最小的那个心磷脂
    }
    #如果所有心磷脂的打分不全为1，则删除分数最小的心磷脂
  }
  #删除错误心磷脂
}

#' @title Scoring a group of Cardiolipin.
#' @param dataframe dataframe(1) including chain,Δ,rt.Its (chain,Δ) can not repeat beside itself.
#' @return dataframe(1) including chain,Δ,rt and their score
#' @seealso [findmyCL::score_main()]
#' @export
score <- function(dataframe) {
  dataframe <-
    apply(X = dataframe , MARGIN = 1 , FUN = as.numeric) %>%
    t() %>%
    as.data.frame()
  colnames(dataframe) <- c("chain" , "Δ" , "rt")
  #将字符型转换为数值型
  data <- findmyCL::turnDataframeList(dataframe = dataframe)
  #转换为列表进行计算
  score_result <-
    purrr::map(.x = data ,
               .f = score_main ,
               dataframe = dataframe) %>%
    plyr::rbind.fill.matrix()
  #进行打分
  colnames(score_result) <- "score"
  #重命名列名
  result <- cbind(dataframe , score_result)
  #列整合结果
  return(result)
}

#' @title Scoring a single Cardiolipin.Its (chain,Δ) can not repeat beside itself.
#' @param x dataframe(1) the single Cardiolipin
#' @param dataframe dataframe(1) including all Cardiolipin
#' @return numeric(1) the score
#' @export
score_main <- function(x , dataframe) {
  if (class(x) == "numeric") {
    x <- x %>%
      t() %>%
      as.data.frame()
  }
  score <- 1
  #分数初始化
  count <- 0
  #记录和多少个其他心磷脂比较过
  same_chain <- which(x$chain == dataframe$chain)
  same_Unsaturation <- which(x$Δ == dataframe$Δ)
  same_chain <- dataframe[same_chain,]
  same_Unsaturation <- dataframe[same_Unsaturation,]
  #相同的链长或者是不饱和度
  for (i in 1:length(same_chain[, 1])) {
    if (same_chain[i,]$Δ  > x$Δ) {
      if (same_chain[i,]$rt < x$rt) {
        score <- score + 1
        count <- count + 1
        #如果符合碳数规律，分数 + 1
      }
      if (same_chain[i,]$rt > x$rt) {
        score <- score - 1
        count <- count + 1
        #如果不符合碳数规律，分数 - 1
      }
    }
    #如果同链长中的不饱和度比 被比较的样本大

    if (same_chain[i,]$Δ  < x$Δ) {
      if (same_chain[i,]$rt > x$rt) {
        score <- score + 1
        count <- count + 1
        #如果符合碳数规律，分数 + 1
      }
      if (same_chain[i,]$rt < x$rt) {
        score <- score - 1
        count <- count + 1
        #如果不符合碳数规律，分数 - 1
      }
    }
    #如果同链长中的不饱和度比 被比较的样本小
  }
  #对链长相同的心磷脂进行计分

  for (i in 1:length(same_Unsaturation[, 1])) {
    if (same_Unsaturation[i,]$chain > x$chain) {
      if (same_Unsaturation[i,]$rt > x$rt) {
        score <- score + 1
        count <- count + 1
        #如果符合碳数规律，分数 + 1
      }
      if (same_Unsaturation[i,]$rt < x$rt) {
        score <- score - 1
        count <- count + 1
        #如果不符合碳数规律，分数 - 1
      }
    }
    #如果同不饱和度中的链长比 被比较的样本大

    if (same_Unsaturation[i,]$chain < x$chain) {
      if (same_Unsaturation[i,]$rt < x$rt) {
        score <- score + 1
        count <- count + 1
        #如果符合碳数规律，分数 + 1
      }
      if (same_Unsaturation[i,]$rt > x$rt) {
        score <- score - 1
        count <- count + 1
        #如果不符合碳数规律，分数 - 1
      }
    }
    #如果同不饱和度中的链长比 被比较的样本小
  }
  #对不饱和度相同的心磷脂进行计分

  score <- score / (count + 1) %>%
    return()
  #最终得分
}
