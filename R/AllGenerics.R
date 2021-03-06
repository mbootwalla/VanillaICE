setGeneric("calculateEmission", function(object, ...) standardGeneric("calculateEmission"))
setGeneric("computeEmission", function(object, hmmOptions) standardGeneric("computeEmission"))
setGeneric("computeHmm", function(object, hmmOptions) standardGeneric("computeHmm"))
setGeneric("emissionPr", function(object) standardGeneric("emissionPr"))
setGeneric("emissionPr<-", function(object, value) standardGeneric("emissionPr<-"))
setGeneric("rangedData", function(object) standardGeneric("rangedData"))
setGeneric("rangedData<-", function(object, value) standardGeneric("rangedData<-"))
setGeneric("segmentData", function(object) standardGeneric("segmentData"))
setGeneric("segmentData<-", function(object, value) standardGeneric("segmentData<-"))

##RangedDataCn generics
setGeneric("LLR", function(object) standardGeneric("LLR"))
setGeneric("nMarkers", function(object) standardGeneric("nMarkers"))
setGeneric("state", function(object) standardGeneric("state"))
setGeneric("rd2df", function(object, hmm.params, ...) standardGeneric("rd2df"))
##setGeneric("plot", useAsDefault=function(x,y,...) graphics::plot(x,y,...))
##plot <- function(object, hmm.params, ...) UseMethod("plot")
##plot.default <- graphics:::plot
setGeneric("plot", function(object, hmm.params, ...) standardGeneric("plot"))

setGeneric("hmm", function(object, hmm.params, ...) standardGeneric("hmm"))
setGeneric("hmm2", function(object, hmm.params, ...) standardGeneric("hmm2"))

setGeneric("snpsetClass", function(object) standardGeneric("snpsetClass"))
setGeneric("copynumberStates", function(object) standardGeneric("copynumberStates"))
setGeneric("states", function(object) standardGeneric("states"))
setGeneric("ICE", function(object) standardGeneric("ICE"))
setGeneric("is.log", function(object) standardGeneric("is.log"))
setGeneric("scaleSds", function(object) standardGeneric("scaleSds"))
setGeneric("log.initialPr", function(object) standardGeneric("log.initialPr"))
setGeneric("normalIndex", function(object) standardGeneric("normalIndex"))
setGeneric("prGtHom", function(object) standardGeneric("prGtHom"))
setGeneric("prGtMis", function(object) standardGeneric("prGtMis"))
setGeneric("prHetCalledHom", function(object) standardGeneric("prHetCalledHom"))
setGeneric("prHetCalledHet", function(object) standardGeneric("prHetCalledHet"))
setGeneric("prHomInNormal", function(object) standardGeneric("prHomInNormal"))
setGeneric("prHomInRoh", function(object) standardGeneric("prHomInRoh"))
setGeneric("rohStates", function(object) standardGeneric("rohStates"))
setGeneric("verbose", function(object) standardGeneric("verbose"))
setGeneric("emit", function(object, hmm.params) standardGeneric("emit"))
setGeneric("markerIndex", function(object) standardGeneric("markerIndex"))
setGeneric("sampleIndex", function(object) standardGeneric("sampleIndex"))

setGeneric("tau", function(object) standardGeneric("tau"))
setGeneric("n2a", function(object) standardGeneric("n2a"))
setGeneric("a2n", function(object) standardGeneric("a2n"))
setGeneric("a2a", function(object) standardGeneric("a2a"))


