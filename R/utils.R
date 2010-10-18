getCrlmmReference <- function(x) paste(annotation(x), "Conf", sep="")
isLoaded <- function(dataset, environ=.vanillaIcePkgEnv)
	exists(dataset, envir=environ)
getVarInEnv <- function(dataset, environ=.vanillaIcePkgEnv){
	if (!isLoaded(dataset))
		stop("Variable ", dataset, " not found in .vanillaIcePkgEnv")
	environ[[dataset]]
}

loader <- function(theFile, envir, pkgname){
	theFile <- file.path(system.file(package=pkgname),
			     "extdata", theFile)
	if (!file.exists(theFile))
		stop("File ", theFile, "does not exist in ", pkgname)
	load(theFile, envir=envir)
}
