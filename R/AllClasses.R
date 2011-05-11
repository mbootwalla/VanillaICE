setClass("RangedDataCn", contains="RangedData")
setValidity("RangedDataCn", function(object){
	all(c("chrom", "sampleId", "numMarkers", "LLR") %in% colnames(object))
})

setClass("HmmOptionList",
	 representation(snpsetClass="character",
			copynumberStates="numeric",
			states="character",
			ICE="logical",
			copyNumber="logical",
			is.log="logical",
			scaleSds="logical",
			log.initialPr="numeric",
			normalIndex="integer",
			prGtHom="numeric",
			prGtMis="numeric",
			prHetCalledHom="numeric",
			prHetCalledHet="numeric",
			prHomInNormal="numeric",
			prHomInRoh="numeric",
			rohStates="logical",
			verbose="logical"))



setMethod("initialize", "HmmOptionList",
	  function(.Object, ...){
		  .Object <- callNextMethod()
	  })
setClass("data.frame.CN", contains="data.frame")

##setClass("track", representation(x="numeric", y="numeric"))

