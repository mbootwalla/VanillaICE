setMethod("hmm", signature(object="oligoSnpSet", hmm.params="HmmOptionList"),
	  function(object, hmm.params, ...){
		  ##browser()
		  ##viterbi(object, hmm.params, ...)
		  viterbi(object, hmm.params, ...)
})
setMethod("hmm", signature(object="CopyNumberSet", hmm.params="HmmOptionList"),
	  function(object, hmm.params, ...){
		  ##viterbi(object, hmm.params, ...)
		  viterbi(object, hmm.params, ...)
})
setMethod("hmm", signature(object="SnpSet", hmm.params="HmmOptionList"),
	  function(object, hmm.params, ...){
		  ##viterbi(object, hmm.params, ...)
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
	  function(object, hmm.params, verbose=1, ...){
		  res <- vector("list", ncol(object))
		  v2 <- verbose==2
		  chromosomes <- unique(chromosome(object))
		  ix <- order(chromosome(object), position(object))
		  if(any(diff(ix) < 0)) {
			  object <- object[ix, ]
		  }
		  marker.index.list <- split(seq(length=nrow(object)), chromosome(object))
		  for(j in 1:ncol(object)){
			  tmp <- vector("list", length(chromosomes))
			  if(verbose > 0){
				  if(j %% 25==0) cat("sample", j, "of", ncol(oligoSet), "\n")
			  }
			  for(k in seq_along(chromosomes)){
				  CHR <- chromosomes[k]
				  i <- marker.index.list[[k]]
				  obj <- object[i, j]
				  tmp[[k]] <- hmm(object=obj, verbose=v2, ...)
			  }
			  if(length(tmp) > 1){
				  rdlist <- RangedDataList(tmp)
				  rd <- stack(rdlist)
				  ix <- match("sample", colnames(rd))
				  if(length(ix) > 0) rd <- rd[, -ix]
				  res[[j]] <- rd
				  rm(rd, rdlist); gc()
			  }
		  }
		  rdlist <- RangedDataList(res)
		  rd <- stack(rdlist)
		  ix <- match("sample", colnames(rd))
		  if(length(ix) > 0) rd <- rd[, -ix]
		  res <- stack(rd)
		  return(res)
	  })
