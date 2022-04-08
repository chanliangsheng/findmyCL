#' @title A function to define how the data store in the findmyCL object
#' @export
defineMyClass <- function() {
  class_cache <- new.env(parent = emptyenv())
  #生成一个全局的环境，用于存放设置的类
  setClass(
    Class = "findmyCLclass",
    slots = list(
      path = "character",
      file = "character",
      MS2 = "MSnExp",
      xcms = "xcmsSet",
      loadCentroidData_parameter = "list",
      ms1MatchResult = "list",
      ms2CheckResult = "list",
      noMS2 = "list",
      ms2MatchResult = "list"
    ),
    where = class_cache
  )
  #where指的是环境的位置
}
