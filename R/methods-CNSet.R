setMethod("hmm2", signature(object="CNSet", hmm.params="HmmOptionList"),
	  function(object, hmm.params, ...){
		  sample.index <- sampleIndex(hmm.params)
		  marker.index <- markerIndex(hmm.params)
		  v <- verbose(hmm.params)
		  if(is.null(sample.index)){
			  if(v > 0) message("sample.index is NULL.  Using all samples")
			  sample.index <- seq(length=ncol(object))
		  }
		  sample.index.list <- splitIndicesByLength(sample.index, ocSamples())
		  if(is.null(marker.index)){
			  if(v > 0) message("marker.index is NULL.  Using all markers")
			  marker.index <- seq(length=nrow(object))
		  }
		  marker.index.list <- split(marker.index, chromosome(object)[marker.index])
		  oligoClasses:::open(object$SNR)
		  snr <- object$SNR[]
		  if(any(snr < 5)) warning(sum(snr < 5), " samples have SNR < 5.  May want to exclude these samples by specifying sample.index.")
		  oligoClasses:::open(object$SNR)
		  rm(snr)
		  fit <- vector("list", length(marker.index.list)*length(sample.index.list))
		  k <- 1
		  oligoClasses:::open(object)
		  for(i in seq_along(marker.index.list)){
			  ii <- marker.index.list[[i]]
			  for(j in seq_along(sample.index.list)){
				  jj <- sample.index.list[[j]]
				  cnset <- object[ii, jj]
				  cnset <- cnset[order(position(cnset)), ]
				  oligoSet <- as(cnset, "oligoSnpSet")
				  oligoSet <- centerCopyNumber(oligoSet, ...)
				  sds <- robustSds2(copyNumber(oligoSet), ...)
				  cnConfidence(oligoSet) <- 1/sds
				  rm(sds, cnset); gc()
				  fit[[k]] <- hmm(oligoSet, hmm.params, ...)
				  k <- k+1
			  }
		  }
		  oligoClasses:::close(object)
		  if(length(fit) > 1){
			  rd <- RangedDataList(fit)
			  rd2 <- stack(rd)
			  rm(fit,rd); gc()
			  ix <- match("sample", colnames(rd2))
			  if(length(ix) > 0) rd2 <- rd2[, -ix]
		  } else rd2 <- fit[[1]]
		  return(rd2)
	  })
