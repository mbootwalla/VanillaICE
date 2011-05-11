newHmmOptionList <- function(snpsetClass="oligoSnpSet",
			     copynumberStates=0:4,
			     states=paste("state", 1:length(copynumberStates), sep=""),
			     ICE=FALSE,
			     copyNumber=TRUE,
			     is.log=FALSE,
			     scaleSds=TRUE,
			     log.initialPr=log(rep(1/length(states), length(states))),
			     normalIndex=3L,
			     prGtHom=numeric(), ##only used when ICE=FALSE
			     prGtMis=rep(1/length(states), length(states)), ##not applicable when ICE is TRUE (no NA's from crlmm genotypes)
			     prHetCalledHom=0.001, ## ignored unless ICE is TRUE
			     prHetCalledHet=0.995, ## ignored unless ICE is TRUE
			     prHomInNormal=0.8,    ## ignored unless ICE is TRUE
			     prHomInRoh=0.999, ## ignored unless ICE is TRUE
			     rohStates=logical(), ## ignored unless ICE is TRUE
			     verbose=TRUE,
			     ...){
	new("HmmOptionList",
	    snpsetClass=snpsetClass,
	    copynumberStates=copynumberStates,
	    states=states,
	    ICE=ICE,
	    copyNumber=copyNumber,
	    is.log=is.log,
	    scaleSds=scaleSds,
	    log.initialPr=log.initialPr,
	    normalIndex=normalIndex,
	    prGtHom=prGtHom,
	    prGtMis=prGtMis,
	    prHetCalledHom=prHetCalledHom,
	    prHetCalledHet=prHetCalledHet,
	    prHomInNormal=prHomInNormal,
	    prHomInRoh=prHomInRoh,
	    rohStates=rohStates,
	    verbose=verbose)
}

setMethod("snpsetClass", "HmmOptionList", function(object) object@snpsetClass)
setMethod("copynumberStates", "HmmOptionList", function(object) object@copynumberStates)
setMethod("states", "HmmOptionList", function(object) object@states)
setMethod("ICE", "HmmOptionList", function(object) object@ICE)
setMethod("copyNumber", "HmmOptionList", function(object) object@copyNumber)
setMethod("is.log", "HmmOptionList", function(object) object@is.log)
setMethod("scaleSds", "HmmOptionList", function(object) object@scaleSds)
setMethod("log.initialPr", "HmmOptionList", function(object) object@log.initialPr)
setMethod("prGtHom", "HmmOptionList", function(object) object@prGtHom)
setMethod("prGtMis", "HmmOptionList", function(object) object@prGtMis)
setMethod("prHetCalledHom", "HmmOptionList", function(object) object@prHetCalledHom)
setMethod("prHetCalledHet", "HmmOptionList", function(object) object@prHetCalledHet)
setMethod("prHomInNormal", "HmmOptionList", function(object) object@prHomInNormal)
setMethod("prHomInRoh", "HmmOptionList", function(object) object@prHomInRoh)
setMethod("rohStates", "HmmOptionList", function(object) object@rohStates)
setMethod("verbose", "HmmOptionList", function(object) object@verbose)


##setMethod("initialize", signature(.Object="HmmOptionList"),
##	  function(.Object,
##		   snpsetClass,
##		   copynumberStates,
##		   states,
##		   ICE,
##		   copyNumber,
##		   is.log,
##		   scaleSds,
##		   log.initialPr,
##		   normalIndex,
##		   prGtHom,
##		   prGtMis,
##		   prHetCalledHom,
##		   prHetCalledHet,
##		   prHomInNormal,
##		   prHomInRoh,
##		   rohStates,
##		   verbose){
##		  ##browser()
##		  ##tmp <- callNextMethod()
##	if(missing(is.log)) stop("Must specify whether the copy number is on the log scale using the <is.log> argument.")
##	if(!snpsetClass %in% c("SnpSet", "CopyNumberSet", "oligoSnpSet")){
##		stop("class must be one of SnpSet, CopyNumberSet, or oligoSet")
##	}
##	if(snpsetClass == "SnpSet") copyNumber <- FALSE
##	stopifnot(is.numeric(normalIndex))
##	stopifnot(length(prGtHom) != length(states))
##	if(ICE) stop("not implemented")
##	.Object <- callNextMethod(snpsetClass=snpsetClass,
##				  copynumberStates=copynumberStates,
##				  states=states,
##				  ICE=ICE,
##				  copyNumber=copyNumber,
##				  is.log=is.log,
##				  scaleSds=scaleSds,
##				  log.initialPr=log.initialPr,
##				  normalIndex=normalIndex,
##				  prGtHom=prGtHom,
##				  prGtMis=prGtMis,
##				  prHetCalledHom=prHetCalledHom,
##				  prHetCalledHet=prHetCalledHet,
##				  prHomInNormal=prHomInNormal,
##				  prHomInRoh=prHomInRoh,
##				  rohStates=rohStates,
##				  verbose=verbose)
##})
