.onLoad <- function(lib, pkg){
   library.dynam("recipeB", pkg, lib)
}
.onUnload <- function(libpath)
    library.dynam.unload("recipeB", libpath)
