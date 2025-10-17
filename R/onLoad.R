.onLoad <- function(libname, pkgname) {
  ns <- getNamespace(pkgname)
  S4_register(CompartmentCall, ns)
  S4_register(CompartmapCall, ns)
  S4_register(MultiCompartmentCall, ns)
  S7::methods_register()
}

# https://github.com/tidyverse/dbplyr/blob/0df138ea3dab075cfdb78c17f41fa3e1a9a5be6e/R/zzz.R#L13-L21
# Silence R CMD check note:
# ** checking whether the namespace can be loaded with stated dependencies ... NOTE
# Warning in .undefineMethod("initialize", Class, classWhere) :
#   no generic function 'initialize' found
#
# I'm not sure why this is necessary, but I suspect it's due to the use of
# setOldClass onLoad
#' @importFrom methods initialize
NULL
