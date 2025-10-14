.onLoad <- function(libname, pkgname) {
  ns <- getNamespace(pkgname)
  S4_register(CompartmentCall, ns)
  S7::methods_register()
}
