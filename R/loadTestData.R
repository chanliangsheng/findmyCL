loadTestData <- function(){
  findmyCL::defineMyClass()
  load(file = "D:/findmyCL/data/MS2matchresult.rdata")
  #载入数据
  return(MS2matchresult)
}
