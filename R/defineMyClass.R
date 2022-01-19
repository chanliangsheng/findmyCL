#' @title A function to define how the data store in the findmyCL object
#' @export
defineMyClass <- function(){
  setClass(Class = "findmyCLclass",
           slots = list(path = "character",
                        file = "character",
                        MS2 = "MSnExp",
                        xcms = "xcmsSet",
                        loadCentroidData_parameter = "list",
                        ms1MatchResult = "list",
                        ms2CheckResult = "list",
                        noMS2 = "list",
                        ms2MatchResult = "list"),
           package = "findmyCL")
}
