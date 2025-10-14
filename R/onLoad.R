.onLoad <- function(libname, pkgname) {
  ns <- getNamespace(pkgname)
  S4_register(CompartmentCall, ns)
  S4_register(CompartmapCall, ns)
  S7::methods_register()
}
