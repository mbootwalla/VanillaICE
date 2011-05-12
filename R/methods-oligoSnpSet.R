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

setMethod("hmm2", signature(object="CNSet", hmm.params="HmmOptionList"),
	  function(object, hmm.params, sample.index, marker.index, centerCN, DF.SD, ...){
		  if(missing(sample.index)){
			  message("sample.index missing.  Using all samples")
			  sample.index <- seq(length=ncol(object))
		  }
		  sample.index.list <- splitIndicesByLength(sample.index, ocSamples())
		  if(missing(marker.index)){
			  message("marker.index missing.  Using all markers")
			  marker.index <- seq(length=nrow(object))
		  }
		  marker.index.list <- splitIndicesByLength(marker.index, chromosome(object))
		  open(object$SNR)
		  snr <- object$SNR[]
		  if(any(snr < 5)) warning(sum(snr < 5), " samples have SNR < 5.  May want to exclude these samples by specifying sample.index.")
		  close(object$SNR)
		  rm(snr)
		  fit <- vector("list", length(marker.index.list)*length(sample.index.list))
		  k <- 1
		  for(i in seq_along(marker.index.list)){
			  ii <- marker.index.list[[i]]
			  for(j in seq_along(sample.index.list)){
				  jj <- sample.index.list[[j]]
				  cnset <- object[ii, jj]
				  oligoSet <- as(cnset, "oligoSnpSet")
				  if(missing(centerCN) | centerCN==FALSE){
					  copyNumber(oligoSet) <- center(copyNumber(oligoSet), at=2)
				  }
				  sds <- robustSds2(copyNumber(oligoSet), DF.PRIOR=DF.SD)
				  cnConfidence(oligoSet) <- 1/sds
				  rm(sds, cnSet); gc()
				  fit[[k]] <- hmm(oligoSet, hmm.params, ...)
				  k <- k+1
			  }
		  }
		  rd <- RangedDataList(fit)
		  rd2 <- stack(rd)
		  rm(fit,rd); gc()
		  ix <- match("sample", colnames(rd2))
		  if(length(ix) > 0) rd2 <- rd2[, -ix]
		  return(rd2)
	  })

setMethod("hmm2", signature(object="oligoSnpSet", hmm.params="HmmOptionList"),
	  function(object, hmm.params, center, ...){
		  res <- vector("list", ncol(object))
		  v2 <- verbose(hmm.params)
		  chromosomes <- unique(chromosome(object))
		  ix <- order(chromosome(object), position(object))
		  if(any(diff(ix) < 0)) {
			  object <- object[ix, ]
		  }
		  marker.index.list <- split(seq(length=nrow(object)), chromosome(object))
		  for(j in 1:ncol(object)){
			  tmp <- vector("list", length(chromosomes))
			  if(v2 > 0){
				  if(j %% 25==0) cat("sample", j, "of", ncol(oligoSet), "\n")
			  }
			  for(k in seq_along(chromosomes)){
				  CHR <- chromosomes[k]
				  i <- marker.index.list[[k]]
				  obj <- object[i, j]
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
