.onLoad <- function(libname, pkgname)
{
  rfamClanDefinitions <<- rfamGetClanDefinitions()
}