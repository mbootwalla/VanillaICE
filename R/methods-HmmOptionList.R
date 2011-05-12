newHmmOptionList <- function(object,
			     snpsetClass=class(object),
			     copynumberStates=0:4,
			     states=paste("state", 1:length(copynumberStates), sep=""),
			     ICE=FALSE,
			     copyNumber=TRUE,
			     is.log=FALSE,
			     scaleSds=TRUE,
			     log.initialPr=log(rep(1/length(states), length(states))),
			     normalIndex=3L,
			     prGtHom=c(1/3, 0.99, 0.7, 0.6, 0.6), ##only used when ICE=FALSE
			     prGtMis=rep(1/length(states), length(states)), ##not applicable when ICE is TRUE (no NA's from crlmm genotypes)
			     prHetCalledHom=0.001, ## ignored unless ICE is TRUE
			     prHetCalledHet=0.995, ## ignored unless ICE is TRUE
			     prHomInNormal=0.8,    ## ignored unless ICE is TRUE
			     prHomInRoh=0.999, ## ignored unless ICE is TRUE
			     rohStates=logical(), ## ignored unless ICE is TRUE
			     tau=1e8,
			     a2n=1,
			     n2a=1,
			     a2a=1,
			     marker.index=NULL,
			     sample.index=NULL,
			     verbose=2L,
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
	    tau=tau,
	    a2n=a2n,
	    n2a=n2a,
	    a2a=a2a,
	    marker.index=marker.index,
	    sample.index=sample.index,
	    verbose=verbose)
}
setValidity("HmmOptionList", function(object){
	ice <- ICE(object)
	if(!ice){
		S <- length(states(object))
		check <- S == length(prGtHom(object))
		if(!check) return(FALSE)
		check <- S == length(prGtMis(object))
		if(!check) return(FALSE)
		check <- S == length(copynumberStates(object))
		if(!check) return(FALSE)
		check <- S == length(log.initialPr(object))
		if(!check) return(FALSE)
		check <- sum(exp(log.initialPr(object))) == 1
		if(!check) return(FALSE)
		if(!is.null(markerIndex(object))){
			check <- all(markerIndex(object) %in% seq(length=nrow(object)))
			if(!check) return(FALSE)
		}
		if(!is.null(sampleIndex(object))){
			check <- all(sampleIndex(object) %in% seq(length=ncol(object)))
			if(!check) return(FALSE)
		}
	} else{
		FALSE
	}
	TRUE
})



setMethod("snpsetClass", "HmmOptionList", function(object) object@snpsetClass)
setMethod("copynumberStates", "HmmOptionList", function(object) object@copynumberStates)
setMethod("states", "HmmOptionList", function(object) object@states)
setMethod("ICE", "HmmOptionList", function(object) object@ICE)
setMethod("copyNumber", "HmmOptionList", function(object) object@copyNumber)
setMethod("is.log", "HmmOptionList", function(object) object@is.log)
setMethod("scaleSds", "HmmOptionList", function(object) object@scaleSds)
setMethod("normalIndex", "HmmOptionList", function(object) object@normalIndex)
setMethod("log.initialPr", "HmmOptionList", function(object) object@log.initialPr)
setMethod("prGtHom", "HmmOptionList", function(object) object@prGtHom)
setMethod("prGtMis", "HmmOptionList", function(object) object@prGtMis)
setMethod("prHetCalledHom", "HmmOptionList", function(object) object@prHetCalledHom)
setMethod("prHetCalledHet", "HmmOptionList", function(object) object@prHetCalledHet)
setMethod("prHomInNormal", "HmmOptionList", function(object) object@prHomInNormal)
setMethod("prHomInRoh", "HmmOptionList", function(object) object@prHomInRoh)
setMethod("rohStates", "HmmOptionList", function(object) object@rohStates)
setMethod("verbose", "HmmOptionList", function(object) object@verbose)
setMethod("tau", "HmmOptionList", function(object) object@tau)
setMethod("a2n", "HmmOptionList", function(object) object@a2n)
setMethod("n2a", "HmmOptionList", function(object) object@n2a)
setMethod("a2a", "HmmOptionList", function(object) object@a2a)
setMethod("markerIndex", "HmmOptionList", function(object) object@marker.index)
setMethod("sampleIndex", "HmmOptionList", function(object) object@sample.index)


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

