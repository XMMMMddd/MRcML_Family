.onLoad <- function(libname, pkgname) {
  Rcpp::loadModule("MRcML_family_MA", TRUE)
}
