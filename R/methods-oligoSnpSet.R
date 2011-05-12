setMethod("hmm", signature(object="oligoSnpSet", hmm.params="HmmOptionList"),
	  function(object, hmm.params, ...){
		  log.emission <- emit(object, hmm.params)
		  dimnames(log.emission) <- list(featureNames(object),
						 sampleNames(object),
						 states(hmm.params))
		  viterbi(object, hmm.params, log.E=log.emission)
})
setMethod("hmm", signature(object="CopyNumberSet", hmm.params="HmmOptionList"),
	  function(object, hmm.params, ...){
		  log.emission <- emit(object, hmm.params)
		  viterbi(object, hmm.params, ...)
})
setMethod("hmm", signature(object="SnpSet", hmm.params="HmmOptionList"),
	  function(object, hmm.params, ...){
		  log.emission <- emit(object, hmm.params)
		  viterbi(object, hmm.params, ...)
})


##setMethod("hmm", signature(object="oligoSnpSet", hmm.params="missing"),
##	  function(object, hmm.params, ...){
##		  ##browser()
##		  ##viterbi(object, hmm.params, ...)
##		  hmmoptions <- hmm.setup(object=object, ...)
##		  callGeneric(object, hmmoptions, ...)
##		  ##hmm(object, hmmoptions)
##})
##setMethod("hmm", signature(object="CopyNumberSet", hmm.params="missing"),
##	  function(object, hmm.params, ...){
##		  ##viterbi(object, hmm.params, ...)
##		  hmmoptions <- hmm.setup(object=object, ...)
##		  ##hmm(object, hmmoptions)
##		  callGeneric(object, hmmoptions, ...)
##})
##setMethod("hmm", signature(object="SnpSet", hmm.params="missing"),
##	  function(object, hmm.params, ...){
##		  ##viterbi(object, hmm.params, ...)
##		  hmmoptions <- hmm.setup(object=object, ...)
##		  ##hmm(object, hmmoptions)
##		  callGeneric(object, hmmoptions, ...)
##})



setMethod("hmm2", signature(object="oligoSnpSet", hmm.params="HmmOptionList"),
	  function(object, hmm.params, ...){
		  res <- vector("list", ncol(object))
		  v2 <- verbose(hmm.params)
		  if(is.null(markerIndex(hmm.params))){
			  marker.index.list <- split(seq(length=nrow(object)), chromosome(object))
			  chromosomes <- unique(chromosome(object))
			  ix <- order(chromosome(object), position(object))
		  } else {
			  marker.index <- markerIndex(hmm.params)
			  marker.index.list <- split(marker.index, chromosome(object)[marker.index])
			  chromosomes <- unique(chromosome(object)[marker.index])
			  ix <- order(chromosome(object)[marker.index], position(object)[marker.index])
		  }
		  if(any(diff(ix) < 0)) {
			  object <- object[ix, ]
		  }
		  if(is.null(sampleIndex(hmm.params))){
			  sample.index <- seq(length=ncol(object))
		  } else sample.index <- sampleIndex(object)
		  for(j in seq_along(sample.index)){
			  jj <- sample.index[j]
			  tmp <- vector("list", length(chromosomes))
			  if(v2 > 0){
				  if(j %% 25==0) cat("sample", j, "of", length(sample.index), "\n")
			  }
			  for(k in seq_along(chromosomes)){
				  CHR <- chromosomes[k]
				  i <- marker.index.list[[k]]
				  obj <- object[i, jj]
				  tmp[[k]] <- hmm(object=obj, hmm.params=hmm.params, ...)
			  }
			  if(length(tmp) > 1){
				  rdlist <- RangedDataList(tmp)
				  rd <- stack(rdlist)
				  ix <- match("sample", colnames(rd))
				  if(length(ix) > 0) rd <- rd[, -ix]
				  rm(rdlist)
			  } else rd <- tmp[[1]]
			  res[[j]] <- rd
			  rm(tmp, rd); gc()
		  }
		  if(length(res) > 1){
			  ##rdlist <- lapply(res, function(x) as(x, "RangedData"))
			  rdlist <- RangedDataList(res)
			  rd <- stack(rdlist)
			  ix <- match("sample", colnames(rd))
			  if(length(ix) > 0) rd <- rd[, -ix]
			  rm(rdlist)
		  } else rd <- res[[1]]
		  return(rd)
	  })

setMethod("emit", signature(object="oligoSnpSet", hmm.params="HmmOptionList"),
	  function(object, hmm.params){
		  ICE <- ICE(hmm.params)
		  ##EMIT.THR <- hmm.params[["EMIT.THR"]]
		  states <- states(hmm.params)
		  verbose <- verbose(hmm.params)
		  normalIndex <- normalIndex(hmm.params)
		  if(all(is.na(cnConfidence(object)))){
			  message("cnConfidence missing.  Using MAD")
			  sds <- robustSds2(copyNumber(object))
			  cnConfidence(object) <- 1/sds
		  }
		  log.cn.emission <- calculateEmission.copynumber2(object,
								   hmm.params)
		  if(!ICE){
			  log.gt.emission <- calculateEmission.genotype(object, hmm.params)
		  } else {
			  ##assumed order
			  ## ROH, normal
			  stop('need to update')
			  log.gt.emission <- array(NA, dim=c(nrow(object), ncol(object), length(states)),
						   dimnames=list(featureNames(object),
						   sampleNames(object),
						   states))
			  tmp <- genotypeEmissionCrlmm(object, hmm.params)
			  rohStates <- which(hmm.params[["rohStates"]])
			  notRohState <- which(!hmm.params[["rohStates"]])
			  for(j in rohStates){
				  log.gt.emission[, , j] <- tmp[, , "ROH"]
			  }
			  for(j in notRohState){
				  log.gt.emission[, , j] <- tmp[, , "normal"]
			  }
		  }
		  log.emission <- log.gt.emission+log.cn.emission
		  if(any(is.na(log.emission))){
			  if(verbose==2) message("Converting missing values in the emission matrix to 0")
			  log.emission[is.na(log.emission)] <- 0
		  }
		  return(log.emission)
})
