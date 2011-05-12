setClass("RangedDataCn", contains="RangedData")
setValidity("RangedDataCn", function(object){
	all(c("chrom", "sampleId", "numMarkers", "LLR") %in% colnames(object))
})

setClassUnion("NullOrNumeric", c("NULL", "numeric"))

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
			n2a="numeric",
			a2n="numeric",
			a2a="numeric",
			tau="numeric",
			marker.index="NullOrNumeric",
			sample.index="NullOrNumeric",
			verbose="integer"))



setMethod("initialize", "HmmOptionList",
	  function(.Object, ...){
		  .Object <- callNextMethod()
	  })

setClass("DataFrameCN", representation(row.names="character",
				       names="character"),
	 contains="list")
setMethod("dimnames", "DataFrameCN", function(x) list(x@row.names, names(x)))
setMethod("dim", "DataFrameCN", function(x) c(length(x@row.names), length(x)))



##setClass("track", representation(x="numeric", y="numeric"))

