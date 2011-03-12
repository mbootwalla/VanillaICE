rowMAD <- function(x, y, ...){
	notna <- !is.na(x)
	sds <- 1.4826*rowMedians(abs(x-rowMedians(x, ...)), ...)
	return(sds)
}

robustSds <- function(x, takeLog=FALSE, ...){
	if(!is.matrix(x)) stop("x is not a matrix")
	if(takeLog) x <- log2(x)
	if(ncol(x) > 3){
		sds1 <- rowMAD(x, ...)
		sds1 <- matrix(sds1, nrow(x), ncol(x))
		sds2 <- apply(x, 2, "mad", ...)
		sds2 <- sds2/median(sds2, na.rm=TRUE)
		sds <- t(t(sds1) * sds2)
	} else {
		sds <- apply(x, 2, "mad", ...)
		sds <- matrix(sds, nrow(x), ncol(x), byrow=TRUE)
	}
	return(sds)
}

viterbi.wrapper <- function(log.emission,
			    log.initial,
			    transitionPr,
			    arm,
			    S,
			    T,
			    result,
			    delta,
			    normal2altered,
			    altered2normal,
			    altered2altered,
			    normalIndex,
			    pAA){
	tmp <- list(as.matrix(as.double(as.matrix(log.emission))),##beta
		    as.double(as.matrix(log.initial)),##initialP
		    as.matrix(as.double(transitionPr)),
		    as.integer(arm),##arm
		    as.integer(S),##number of states
		    as.integer(T),##number of loci
		    result,##placeholder for results
		    as.matrix(as.double(delta)),##delta
		    normal2altered=normal2altered,##c1
		    altered2normal=altered2normal,##c2
		    altered2altered=altered2altered,##c3
		    as.integer(normalIndex),
		    as.double(pAA))##normalState
	.C("viterbi",
	   log.emission=tmp[[1]],
	   log.initial=tmp[[2]],
	   tau=tmp[[3]],
	   arm=tmp[[4]],
	   S=tmp[[5]],
	   T=tmp[[6]],
	   viterbiSeq=tmp[[7]],
	   delta=tmp[[8]],
	   N2A=tmp[[9]],
	   A2N=tmp[[10]],
	   A2A=tmp[[11]],
	   normalIndex=tmp[[12]],
	   pAA=tmp[[13]])  ##can verify that the transition prob. matrix is correct (for last snp)
}

getChromosomeArm <- function(object){
	chrom <- chromosome(object)
	pos <- position(object)
	if(!is.integer(chrom)) {
		chrom <- chromosome2integer(chrom)
	}
	if(!all(chrom %in% 1:24)){
			warning("Chromosome annotation is currently available for chromosomes 1-22, X and Y")
			message("Please add/modify data(chromosomeAnnotation, package='SNPchip') to accomodate special chromosomes")
			stop()
	}
	if(!is.integer(pos)) {
		message("Coerced pos to an integer.")
		pos <- as.integer(pos)
	}
	data(chromosomeAnnotation, package="SNPchip", envir=environment())
	chromosomeAnnotation <- as.matrix(chromosomeAnnotation)
	chrAnn <- chromosomeAnnotation
	uchrom <- unique(SNPchip:::integer2chromosome(chrom))
	chromosomeArm <- vector("list", length(uchrom))
	positionList <- split(pos, chrom)
	positionList <- positionList[match(unique(chrom), names(positionList))]
	for(i in seq(along=uchrom)){
		chromosomeArm[[i]] <- as.integer(ifelse(positionList[[i]] <= chrAnn[uchrom[i], "centromereEnd"], 0, 1))
	}
	chromosomeArm <- unlist(chromosomeArm)
	chromosomeArm <- cumsum(c(0, diff(chromosomeArm) != 0 | diff(chrom) != 0))
	chromosomeArm <- chromosomeArm+1 ##start at 1
}

viterbi <- function(object,
		    hmm.params,
		    verbose=TRUE,
		    normal2altered=1,
		    altered2normal=1,
		    altered2altered=1,
		    TAUP=1e8){
	log.E <- hmm.params[["log.emission"]]
	sns <- colnames(log.E)
	if(is.null(sns)) stop("no dimnames for log.emission")
	log.initial <- hmm.params[["log.initial"]]

	if(normal2altered <= 0) stop("normal2altered must be > 0")
	if(altered2normal <= 0) stop("altered2normal must be > 0")
	if(altered2altered <= 0) stop("altered2altered must be > 0")

	arm <- getChromosomeArm(object)
	normalIndex <- hmm.params[["normalIndex"]]
	if(normalIndex < 1 | normalIndex > dim(log.E)[[3]]){
		stop("normalIndex in hmm.params not valid")
	}
	c1 <- normal2altered
	c2 <- altered2normal
	c3 <- altered2altered
	TT <- nrow(object)
	states <- hmm.params[["states"]]
	S <- length(states)
	delta <- matrix(as.double(0), nrow=TT, ncol=S)
	rangedData <- list()
	for(j in 1:ncol(log.E)){
		rD <- vector("list", length(unique(arm)))
		for(a in seq(along=unique(arm))){
			I <- arm == a
			if(sum(I) < 2) next()
			T <- sum(I)
			transitionPr <- exp(-2 * diff(position(object)[I])/TAUP)
			##is the lower bound a function of normal2altered, altered2normal, altered2altered?
			minimum <- 1-1/((S-1)*c1) + 0.01
			transitionPr[transitionPr < minimum] <- minimum
			if(any(transitionPr < 0 | transitionPr > 1)) stop("Transition probabilities not in [0,1].  Order object by chromosome and physical position")
			result <- rep(as.integer(0), T)
			viterbiResults <- viterbi.wrapper(log.emission=log.E[I, j, ],
							  log.initial=log.initial,
							  transitionPr=transitionPr,
							  arm=arm[I],
							  S=S,
							  T=T,
							  result=result,
							  delta=delta[I, ],
							  normal2altered=normal2altered,
							  altered2normal=altered2normal,
							  altered2altered=altered2altered,
							  normalIndex=normalIndex,
							  pAA=rep(0, S^2))
			M <- matrix(viterbiResults[["pAA"]], S, S)
			if(!all(is.finite(M))) stop("Infinite values in transition prob. matrix. Check that rows are ordered by physical position")
			if(!all.equal(rowSums(exp(M)), rep(1, S))){
				warning("Rows of the transition probability matrix do not sum to 1")
			}
			viterbiSequence <- viterbiResults[["viterbiSeq"]]
			rl <- Rle(viterbiSequence)
			starts <- start(rl)
			LLR <- rep(999,  length(starts))
			log.emission <- matrix(viterbiResults[["log.emission"]], T, S)
			##** The NA is to stagger the transition probabilities by 1
			##  -- this way, the same index can be used to multiply the transition and emission probabilities
			p <- c(NA, as.numeric(viterbiResults[["tau"]]))
			lP.N2N <- log(1-((1-p)*(S-1)*c1)) ##probability normal -> normal
			lP.N2A <- log((1-p)*c1) ##probability normal -> altered
			P.A2A <- sapply(1-((1-p)*(c2+(S-2)*c3)), function(x) max(x, 0.01))
			lP.A2A <- log(P.A2A) ## probability altered to same altered state
			lP.A2N <- log((1-p)*c2) ##probability altered -> normal
			lP.A2Astar <- log((1-p)*c3) ## probability altered -> different altered state
			##For each seqment, compute the likelihood ratio
			names(log.initial) <- hmm.params[["states"]]
			for(k in seq(along=starts)){
				##if(any(LLR < 0)) browser()
				index <- start(rl)[k]:end(rl)[k]
				thisState <- unique(viterbiSequence[index])
				if(thisState == normalIndex){
					LLR[k] <- 0
					next()
				}
				first.index <- min(index)
				last.index <- max(index)
				## 6 Rules  (1 < t < T)
				## 1.  index=1
				## 2.  index=t
				## 3.  index=t,t+1
				## 4.  index=T
				## 5.  index=1,2
				## 6,  index=T-1, T
				##------
				## 1. index=1
				if(first.index == last.index & last.index==1){
					##note the last term cancels
					logLik.vit <- log.initial[thisState]+log.emission[1, thisState] + lP.A2N[2] + log.emission[2, normalIndex]
					logLik.null <- log.initial[normalIndex]+log.emission[1, normalIndex] + lP.N2N[2] + log.emission[2, normalIndex]
					LLR[k] <- logLik.vit-logLik.null
					next()
				}
				##2 index = t
				if(length(index) == 1 & first.index > 1 & last.index < T){
					##note the last term cancels
					logLik.vit <- sum(lP.N2A[index] + log.emission[index, thisState]) + lP.A2N[last.index+1]+log.emission[last.index+1, normalIndex]
					logLik.null <- sum(lP.N2N[index] + log.emission[index, normalIndex]) + lP.N2N[last.index+1] + log.emission[last.index+1, normalIndex]
					LLR[k] <- logLik.vit-logLik.null
					next()
				}
				##if T=t+1?
				## 3: index = t, ..., t*  , t>1, t* < T, t* != t
				if(first.index != last.index & first.index > 1 & last.index < T){
					index2 <- index[-1]
					logLik.vit <- lP.N2A[first.index] +
					              sum(lP.A2A[index2]) +
					              lP.A2N[last.index+1] +
						      log.emission[first.index, thisState] +
						      sum(log.emission[index2, thisState]) +
					              log.emission[last.index+1, normalIndex]
					logLik.null <-
						sum(lP.N2N[index]) +
						lP.N2N[last.index+1]  +
						sum(log.emission[index, normalIndex]) +
						log.emission[last.index+1, normalIndex]
					LLR[k] <- logLik.vit-logLik.null
					next()
				}
				## 4: index = T
				if(first.index == last.index & last.index == T){
					logLik.vit <- lP.N2A[T] + log.emission[T, thisState]
					logLik.null <- lP.N2N[T] + log.emission[T, normalIndex]
					LLR[k] <- logLik.vit-logLik.null
					next()
				}
				## 5: index = 1, 2, ...
				if(first.index != last.index & first.index == 1 & last.index < T){
					index2 <- index[-1]## t=2, ...., t*
					logLik.vit <- log.initial[thisState] + log.emission[first.index, thisState]  + sum(lP.A2A[index2] + log.emission[index2, thisState]) + lP.A2N[last.index+1] + log.emission[last.index+1, normalIndex]
					logLik.null <- log.initial[normalIndex] + log.emission[first.index, normalIndex] + sum(lP.N2N[index2] + log.emission[index2, normalIndex]) + lP.N2N[last.index+1] + log.emission[last.index+1, normalIndex]
					LLR[k] <- logLik.vit-logLik.null
					next()
				}
				if(first.index != last.index & first.index == 1 & last.index == T){
					index2 <- index[-1]## t=2, ...., t*
					logLik.vit <- log.initial[thisState] + log.emission[first.index, thisState]  + sum(lP.A2A[index2] + log.emission[index2, thisState]) ##+ lP.A2N[last.index+1] + log.emission[last.index+1, normalIndex]
					logLik.null <- log.initial[normalIndex] + log.emission[first.index, normalIndex] + sum(lP.N2N[index2] + log.emission[index2, normalIndex]) ##+ lP.N2N[last.index+1] + log.emission[last.index+1, normalIndex]
					LLR[k] <- logLik.vit-logLik.null
					next()
				}
				## 6: index = t, ...T
				if(first.index != last.index & last.index == T){
					index2 <- index[-1]
					logLik.vit <- lP.N2A[first.index] + log.emission[first.index, thisState] + sum(lP.A2A[index2] + log.emission[index2, thisState])
					logLik.null <- lP.N2N[first.index] + log.emission[first.index, normalIndex] + sum(lP.N2N[index2] + log.emission[index2, normalIndex])
					LLR[k] <- logLik.vit - logLik.null
				}
			}
			start.index <- start(rl)
			end.index <- end(rl)
			pos <- position(object)[I]
			##this is tricky since we've added an index to force a segment for each arm.
			start <- pos[start.index]
			end <- pos[end.index]
			##numMarkers <- unlist(numMarkers)
			numMarkers <- width(rl)
			states <- viterbiSequence[start.index]
			ir <- IRanges(start=start, end=end)
			rD[[a]] <- RangedData(ir,
					      ##space=rep(paste("chr", unique(chromosome(object)[I]), sep=""), length(ir)),
					      chrom=rep(unique(chromosome(object)[I]), length(ir)),
					      ##sampleId=sampleNames(object)[j],
					      sampleId=colnames(log.E)[j],
					      state=states,
					      numMarkers=numMarkers,
					      LLR=LLR)
		}
		notnull <- !sapply(rD, is.null)
		rD <- rD[notnull]
		L <- sapply(rD, nrow)
		if(any(L == 1) & any(L > 1)){
			rangedData[[j]] <- c(do.call(c, rD[L == 1]), do.call(c, rD[L > 1]))
		} else {
			rangedData[[j]] <- do.call(c, rD)
		}
	}
##	L <- sapply(rangedData, nrow)
##	if(any(L == 1) & any(L > 1)){
##		rangedData <- c(do.call(c, rangedData[L == 1]), do.call(c, rangedData[L > 1]))
##	} else {
##		rangedData <- do.call(c, rangedData)
##	}
	sampleId <- unlist(lapply(rangedData, function(x) x$sampleId))
	state <- unlist(lapply(rangedData, function(x) x$state))
	numMarkers <- unlist(lapply(rangedData, function(x) x$numMarkers))
	LLR <- unlist(lapply(rangedData, function(x) x$LLR))
	chr <- unlist(lapply(rangedData, function(x) x$chrom))
	starts <- unlist(lapply(rangedData, start))
	ends <- unlist(lapply(rangedData, end))
	ir <- IRanges(start=starts, end=ends)
	rangedData <- RangedData(ir,
				 chrom=chr,
				 sampleId=sampleId,
				 state=state,
				 numMarkers=numMarkers,
				 LLR=LLR)
	return(rangedData)
}






##viterbi <- function(emission,
##		    tau,
##		    arm,
##		    initialStateProbs,
##		    verbose=FALSE,
##		    chromosome,
##		    position,
##		    sampleNames,
##		    locusNames,
##		    normalIndex,
##		    returnLikelihood=FALSE,
##		    normal2altered=1,
##		    altered2normal=1,
##		    altered2altered=1){
##	if(class(emission) != "array") stop("emission probabilities must be an array: snps, samples, states. ")
##	if(missing(sampleNames)) sampleNames <- colnames(emission)
##	if(missing(normalIndex)) stop("Must specify integer for normalIndex")
##	if(!is.numeric(normalIndex)) stop("normalIndex should be numeric")
##	viterbiSequence <- matrix(NA, nrow(emission), ncol(emission), dimnames)
##	S <- dim(emission)[3]
##	T <- nrow(emission)
##	if(missing(initialStateProbs)){
##		initialStateProbs <- log(rep(1/S, S))
##	}
##	if(length(initialStateProbs) != S){
##		stop("initialStateProbs (the initial state probabilities, should be a numeric vector of length S, where S is the number of hidden states")
##	}
##	if(!all(initialStateProbs <= 0)){
##		if(all(initialStateProbs >= 0 & initialStateProbs <= 1)){
##			initialStateProbs <- log(initialStateProbs)
##		} else stop("initial state probabilities should be a probability or a log probability")
##	}
##	if(any(is.na(emission))){
##		if(verbose) message("Converting missing values in the emission matrix to 0")
##		emission[is.na(emission)] <- 0
##	}
##	if(any(is.nan(emission))){
##		message("some of the log emission probabilities are NaN.  Replacing with 0")
##		emission[is.nan(emission)] <- 0
##	}
##	if(any(emission < -50)){
##		message("some of the log emission probabilities are very small -- probable outliers.  Replacing with a small value (-10)")
##		emission[emission < -50] <- -50
##	}
##	if(missing(arm)){
##		message("chromosome arm not specified...HMM is not fit separately to each chromosomal arm")
##		arm <- rep(as.integer(1), T)
##	}
##	if(length(arm) != T) {
##		message("arm not the right length.  assuming all values on same chromosomal arm")
##		arm <- rep(as.integer(1), T)
##	}
##	if(missing(tau)){
##		stop("transition probabilities not specified")
##	}
##	if(length(tau) != T) stop("tau must have length T")
##	## The last index is arbitrary and, by default, is NA. Must replace this by a number-- C can not handle NA's
##	tau[is.na(tau)] <- 0
##	delta <- matrix(as.double(0), nrow=T, ncol=S)
##	browser()
##	rangedData <- list()
##	for(j in 1:ncol(results)){
##		rD <- vector("list", length(unique(arm)))
##		for(a in seq(along=unique(arm))){
##			I <- arm == a
##			T <- sum(I)
##			result <- rep(as.integer(0), T)
##			tmp <- list(as.matrix(as.double(as.matrix(emission[I, j, ]))),##beta
##				    as.double(as.matrix(initialStateProbs)),##initialP
##				    as.matrix(as.double(tau[I])),##tau
##				    as.integer(arm[I]),##arm
##				    as.integer(S),##number of states
##				    as.integer(T),##number of loci
##				    result,##placeholder for results
##				    as.matrix(as.double(delta[I, ])),##delta
##				    normal2altered=normal2altered,##c1
##				    altered2normal=altered2normal,##c2
##				    altered2altered=altered2altered,##c3
##				    as.integer(normalIndex),
##				    as.double(rep(0, S^2)))##normalIndex
##			tmp2 <- .C("viterbi",
##				   emission=tmp[[1]],
##				   initialStateProbs=tmp[[2]],
##				   tau=tmp[[3]],
##				   arm=tmp[[4]],
##				   S=tmp[[5]],
##				   T=tmp[[6]],
##				   viterbiSeq=tmp[[7]],
##				   delta=tmp[[8]],
##				   N2A=tmp[[9]],
##				   A2N=tmp[[10]],
##				   A2A=tmp[[11]],
##				   normalIndex=tmp[[12]],
##				   pAA=tmp[[13]])  ##can verify that the transition prob. matrix is correct (for last snp)
##			##check transition probs.
##			M <- matrix(tmp2[["pAA"]], S, S)
##			if(!all(is.finite(M))) stop("Infinite values in transition prob. matrix")
##			if(!all.equal(rowSums(exp(M)), rep(1, S))){
##				warning("Rows of the transition probability matrix do not sum to 1")
##			}
##			viterbiSequence <- tmp2[["viterbiSeq"]]
##			logInitialP <- initialStateProbs
##			rl <- Rle(viterbiSequence)
##			starts <- start(rl)
##			LLR <- rep(NA,  length(starts))
##			logE <- matrix(tmp2[["emission"]], T, S)
##			p <- as.numeric(tmp2[["tau"]])
##			c1 <- normal2altered
##			c2 <- altered2normal
##			c3 <- altered2altered
##			lP.N2N <- log(1-((1-p)*(S-1)*c1)) ##probability normal -> normal
##			lP.N2A <- log((1-p)*c1) ##probability normal -> altered
##			lP.A2A <- log(1-((1-p)*(c2+(S-2)*c3))) ## probability altered to same altered state
##			lP.A2N <- log((1-p)*c2) ##probability altered -> normal
##			lP.A2Astar <- log((1-p)*c3) ## probability altered -> different altered state
##			##For each seqment, compute the likelihood ratio
##			for(k in seq(along=starts)){
##				index <- start(rl)[k]:end(rl)[k]
##				thisState <- unique(viterbiSequence[index])
##				first.index <- min(index)
##				last.index <- max(index)
##				if(last.index < T){
##					next.index <- last.index+1
##					next.state <- viterbiSequence[next.index]
##					lE.next.state <- logE[next.index, next.state]
##					tp.last <- lP.A2N[last.index] ##not next.index!
##				} else{ ##max(index) = T
##					lE.next.state <- NULL
##					tp.last <- NULL
##				}
##				if(min(index) > 1){
##					tps <- c(lP.N2A[min(index)-1], lP.A2A[index[-length(index)]], tp.last)
##					tps.normal <- lP.N2N[c(index-1, max(index))]
##				} else {
##					tps <- c(logInitialP[thisState], lP.A2A[index[-length(index)]], tp.last)
##					tps.normal <- c(logInitialP[normalIndex],
##							lP.N2N[index[-length(index)]])
##				}
##				lEs <- c(logE[index, thisState], lE.next.state)
##				lE.normal <- c(logE[index, normalIndex], lE.next.state)
##				loglik.Vit <- sum(lEs+tps)
##				loglik.null <- sum(lE.normal+tps.normal)
##				LLR[k] <- loglik.Vit-loglik.null
##			}
##			start.index <- start(rl)
##			end.index <- end(rl)
##			pos <- position[I]
##			##this is tricky since we've added an index to force a segment for each arm.
##			start <- pos[start.index]
##			end <- pos[end.index]
##			##numMarkers <- unlist(numMarkers)
##			numMarkers <- width(rl)
##			states <- viterbiSequence[start.index]
##			ir <- IRanges(start=start, end=end)
##			rD[[a]] <- RangedData(ir,
##					      space=rep(paste("chr", unique(chromosome[I]), sep=""), length(ir)),
##					      sampleId=sampleNames[j],
##					      state=states,
##					      numMarkers=numMarkers,
##					      LLR=LLR)
##		}
##		tmp <- do.call(c, rD[sapply(rD, nrow) == 1])
##		tmp2 <- do.call(c, rD[sapply(rD, nrow) > 1])
##		rD <- c(tmp, tmp2)
##		rangedData[[j]] <- rD
##		##results <- list(stateSequence=results, logLikelihoodRatio=lrdiff)
##	}
##	rangedData <- do.call(c, rangedData)
##	return(rangedData)
##}
trioOptions <- function(opts,
			states=c("BPI", "notBPI"),
			##initialP=c(0.99, 0.01),
			##TAUP=1e7,
			prGtError=c(0.001, 0.01),
			##verbose=FALSE,
			allowHetParent=FALSE,
			##normalIndex=1,
			##normal2altered=1,
			##altered2normal=1,
			useCrlmmConfidence=FALSE){
	names(prGtError) <- states
	opts[["states"]] <- states
	opts$prGtError <- prGtError
	opts$allowHetParent <- allowHetParent
	##useCrlmmConfidence=useCrlmmConfidence)
	return(opts)
}

##setMethod("coerce", c("CNSet", "list"), function(from, to){
##	if(ncol(object) < 3){
##		return("No complete trios")
##	}
##	stopifnot(all(c("familyId", "fatherId", "motherId", "individualId") %in% varLabels(object)))
##
##	for(i in 1:nrow(trios)){
##
##
##	}
##})

hmm.setup <- function(object,
		      states=paste("state", 1:length(copynumberStates), sep=""),
		      ICE=FALSE,
		      copyNumber=TRUE,
		      copynumberStates=0:4,
		      EMIT.THR=-10,
		      scaleSds=TRUE,
		      verbose=TRUE,
		      log.initial=log(rep(1/length(states), length(states))),
		      normalIndex=3,
		      prGenotypeHomozygous=numeric(), ##only used when ICE=FALSE
		      prGenotypeMissing=rep(1/length(states), length(states)), ##not applicable when ICE is TRUE (no NA's from crlmm genotypes)
		      pHetCalledHom=0.001, ## ignored unless ICE is TRUE
		      pHetCalledHet=0.995, ## ignored unless ICE is TRUE
		      pHomInNormal=0.8,    ## ignored unless ICE is TRUE
		      pHomInRoh=0.999, ## ignored unless ICE is TRUE
		      rohStates=logical(), ## ignored unless ICE is TRUE
		      trioHmm=FALSE,
		       ...){  ## whether the save the emission probabilities
	if(!class(object) %in% c("SnpSet", "CopyNumberSet", "oligoSnpSet")){
		message("class of object must be one of SnpSet, CopyNumberSet, or oligoSet")
		stop()
	}
	if(class(object) == "SnpSet") copyNumber <- FALSE
	stopifnot(is.numeric(normalIndex))
	if(ICE){
		##check with crlmm confidence scores for this platform are available
		if(length(annotation(object)) < 1){
			stop("When ICE is TRUE, annotation slot in the SnpSet object must be specified")
		}
		if(length(rohStates) != length(states)){
			stop("Must specify which states are 'ROH-like.  See documentation.")
		}
	}
	opts <- list(copynumberStates=copynumberStates,
		     states=states,
		     ICE=ICE,
		     copyNumber=copyNumber,
		     EMIT.THR=EMIT.THR,
		     scaleSds=scaleSds,
		     log.initial=log.initial,
		     normalIndex=normalIndex,
		     prGenotypeHomozygous=prGenotypeHomozygous,
		     prGenotypeMissing=prGenotypeMissing,
		     pHetCalledHom=pHetCalledHom,
		     pHetCalledHet=pHetCalledHet,
		     pHomInNormal=pHomInNormal,
		     pHomInRoh=pHomInRoh,
		     rohStates=rohStates,
		     log.emission=NULL,
		     log.gt.emission=NULL,
		     log.cn.emission=NULL,
		     verbose=verbose,
		     ...)
	if(trioHmm){
		opts <- trioOptions(opts)
		offspringId <- sampleNames(object)[object$fatherId != 0 & object$motherId != 0]
		trios <- as.matrix(t(sapply(offspringId, findFatherMother, object=object)))
		trios <- trios[rowSums(is.na(trios)) == 0, , drop=FALSE]
		colnames(trios) <- c("father", "mother", "offspring")
		opts$trios <- trios
		opts$copyNumber <- FALSE
		##trioList <- as(object, "TrioSetList")
	}
	##opts[["tau"]] <- transitionProbability(object, opts)
	if(opts$copyNumber){
		##check that the median copy number is near the
		##specified copy number for the normal hidden state
		##(helps prevent errors from forgetting to take the
		##log)
		isAutosome <- chromosome(object) <= 22
		if(sum(isAutosome) > 1){
			med <- median(as.numeric(copyNumber(object)[isAutosome, ]), na.rm=TRUE)
			med.normal <- copynumberStates[normalIndex]
			delta <- abs(med-med.normal)
			if(delta < 0.1){
				message("The absolute difference between the median copy number and copynumberState[normalIndex] is ", abs(med-med.normal))
			} else {
				warning("The absolute difference between the median copy number and copynumberState[normalIndex] is ", delta)
			}
		}
	}
	if(!trioHmm){
		message("Computing emission probabilities.")
		if(!ICE & is(object, "oligoSnpSet")){
			if(length(prGenotypeHomozygous) != length(states)){
				stop("The probability of an AA or BB genotype must be specified for each state.  This specification is obtained by setting by passing a numeric vector of probabilities to the prGenotypeHomozygous argument.")
			}
		}
		opts <- calculateEmission(object, opts)
	}
	opts
}

##transitionProbability.oligoSnpSet <- function(object, hmmOptions){
##	chrom <- chromosome(object)
##	pos <- position(object)
##	TAUP <- hmmOptions[["TAUP"]]
##	calculateTransitionProbability(chrom, pos, TAUP)
##}
##
##setMethod("transitionProbability", "SnpSet", function(object, hmmOptions){
##	transitionProbability.oligoSnpSet(object, hmmOptions)
##})

##calculateTransitionProbability <- function(chromosome, position, TAUP=1e8, chromosomeAnnotation, verbose=FALSE){
##	if(!is.integer(chromosome)) {
##		chromosome <- chromosome2integer(chromosome)
##	}
##	if(!all(chromosome %in% 1:24)){
##			warning("Chromosome annotation is currently available for chromosomes 1-22, X and Y")
##			message("Please add/modify data(chromosomeAnnotation, package='SNPchip') to accomodate special chromosomes")
##			stop()
##	}
##	if(!is.integer(position)) {
##		if(verbose) message("Coerced position to an integer.")
##		position <- as.integer(position)
##	}
##	if(length(chromosome) != length(position)) stop("chromosome and position arguments must be the same length")
##	if(missing(chromosomeAnnotation)){
##		if(verbose) message("chromosomeAnnotation not specified... using centromere locations from SNPchip")
##		data(chromosomeAnnotation, package="SNPchip", envir=environment())
##		chromosomeAnnotation <- as.matrix(chromosomeAnnotation)
##	}
##	chrAnn <- chromosomeAnnotation
##	uchrom <- unique(integer2chromosome(chromosome))
##	chromosomeArm <- tau <- vector("list", length(uchrom))
##	positionList <- split(position, chromosome)
##	positionList <- positionList[match(unique(chromosome), names(positionList))]
##	for(i in seq(along=uchrom)){
##		chromosomeArm[[i]] <- as.integer(ifelse(positionList[[i]] <= chrAnn[uchrom[i], "centromereEnd"], 0, 1))
##		##probability SNP is informative.  Note, the last index is assigned zero -- this is arbitrary.
##		tau[[i]] <- c(exp(-2 * diff(positionList[[i]])/TAUP), NA)
##		if(any(tau[[i]] > 1, na.rm=TRUE)){
##			##values greater than one occur when the diff is negative.
##			stop("check that the physical position is ordered from smallest to largest")
##		}
##		if(any(tau[[i]] < 0, na.rm=TRUE)) stop("some of the computed transition probabilities less than zero.")
##	}
##	chromosomeArm <- unlist(chromosomeArm)
##	chromosomeArm <- cumsum(c(0, diff(chromosomeArm) != 0 | diff(chromosome) != 0))
##	chromosomeArm <- chromosomeArm+1 ##start at 1
##	tau <- unlist(tau)
##	annotation <- cbind(chromosome, position, chromosomeArm, tau)
##	colnames(annotation) <- c("chromosome", "position", "arm", "transitionPr")
##	range.tau <- range(annotation[, "transitionPr"], na.rm=TRUE)
##	##check range
##	if(range.tau[1] < 0 | range.tau[2] > 1) stop("Transition probabilities are not valid.  Check that the snpset object has been ordered by chromosome and physical position.")
##	##assign a minimum pr. that the snp is informative
##	tau[tau < 0.5] <- 0.5
##	annotation[, "transitionPr"] <- tau
##	return(annotation)
##}


##Code this in C
viterbiR <- function(emission, initialP, tau, arm){
##	emission: log emission probabilities
##	initialP: log initial state probabilities (vector of length S)
##	tau: transition probabilities (original scale)
##	tau.scale: matrix on original scale
	S <- ncol(emission)
	T <- nrow(emission)  ##matrix of T rows and S columns
	delta <- psi <- matrix(NA, T, S)
	delta[1, ] <- initialP + emission[1, ]
	psi[1, ] <- rep(0, S)
	i <- which(rowSums(is.na(emission)) > 0)
	tau.scale <- 1
	#emission[i, ] <- 0
	for(t in 2:T){
		if(t %% 10000 == 0) cat(".")
		if(arm[t] != arm[t-1]){
			delta[t, ] <- initialP + emission[t, ]
			psi[t, ] <- 0
			next()
		}
##		AA <- matrix(NA, nr=S, nc=S)
##		for(i in 1:S){
		AA <- matrix(tau[t-1], nr=S, nc=S)
		epsilon <- (1-tau[t-1])/(S-1)
		##eNormal <- (1-tauNormal[t-1])/(S-1)
		AA[upper.tri(AA)] <- AA[lower.tri(AA)] <- epsilon
		##AA[normalIndex, ] <- rep(eNormal, S)
		##AA[normalIndex, normalIndex] <- tauNormal[t-1]
		AA <- log(AA*tau.scale)
		for(j in 1:S){
			tmp <- delta[t-1, ] + AA[, j]
			delta[t, j] <- max(tmp) + emission[t, j]
			psi[t, j] <- order(tmp, decreasing=TRUE)[1]
		}
	}
	Pstar <- max(delta[nrow(delta), ])
	qhat <- rep(NA, nrow(delta))
	qhat[T] <- order(delta[T, ], decreasing=TRUE)[1]
	for(t in (T-1):1){
		if(arm[t] != arm[t+1]){
			qhat[t] <- order(delta[t, ], decreasing=TRUE)[1]
		} else {
			qhat[t] <- psi[t+1, qhat[t+1]]
		}
	}
	return(qhat)
}

##setMethod("hmm", "oligoSnpSet", function(object, hmmOptions){
##	viterbi(object, hmmOptions)
##})
##
##setMethod("hmm", "CNSet", function(object, hmmOptions){
##	viterbi(object, hmmOptions)
##})
##
##setMethod("hmm", "CopyNumberSet", function(object, hmmOptions){
##	viterbi(object, hmmOptions)
##})
hmm <- function(object, hmm.params, ...){
	if(missing(hmm.params)) stop("missing hmm.params.  See hmm.setup")
	if(missing(object)) stop("missing object.")
	viterbi(object, hmm.params, ...)
}


##hmm <- function(object,
##		states,
##		mu=NULL,
##		probs=NULL,
##		takeLog=FALSE,
##		initialP,
##		returnSegments=TRUE,
##		TAUP=1e8,
##		verbose=FALSE,
##		ice=FALSE,
##		envir,
##		normalIndex){
##	if(missing(envir)) envir <- new.env()
##	if(!all(c("position", "chromosome") %in% fvarLabels(object))){
##		stop("'position' and 'chromosome' must be in fvarLabels(object), or transitionPr must be provided")
##	}
##	object <- object[order(chromosome(object), position(object)), ]
##	envir[["locusset"]] <- object
##	if(missing(states)) stop("must specify states")
##	envir[["states"]] <- states
##	if(missing(initialP)) initialP <- rep(1, length(states))/length(states)
##	envir[["initialP"]] <- initialP
##	envir[["mu"]] <- mu
##	if(is.null(probs)){
##		if(ice){
##			probs <- c(0.05, 0.99, 0.7, 0.999)
##		}
##	}
##	envir[["probs"]] <- probs
##	envir[["takeLog"]] <- takeLog
##	envir[["returnSegments"]] <- returnSegments
##	envir[["TAUP"]] <- TAUP
##	envir[["verbose"]] <- verbose
##	envir[["ice"]] <- ice
##	tau <- transitionProbability(chromosome=chromosome(object),
##				     position=position(object),
##				     TAUP=TAUP,
##				     verbose=verbose)
##	arm <- tau[, "arm"]
##	transitionPr <- tau[, "transitionPr"]
##	envir[["transitionPr"]] <- transitionPr
##	envir[["arm"]] <- arm
##	if(takeLog){
##		copyNumber(object) <- log2(copyNumber(object))
##		mu <- log2(mu)
##	}
##	if(verbose) message("Calculating emission probabilities")
##	calculateEmission(object=object,
##			  mu=mu,
##			  probs=probs,
##			  envir=envir,
##			  states=states,
##			  verbose=verbose,
##			  ice=ice)
##	if(!is.null(envir[["emission.gt"]])){
##		emission <- envir[["emission.cn"]]+envir[["emission.gt"]]
##	} else {
##		emission <- envir[["emission.cn"]]
##	}
##	envir[["emission"]] <- emission
##	##emission <- .GlobalEnv[["emission.cn"]]+.GlobalEnv[["emission.gt"]]
##	viterbiResults <- viterbi(initialStateProbs=log(initialP),
##				  emission=emission,
##				  tau=transitionPr,
##				  arm=arm,
##				  verbose=verbose,
##				  normalIndex=normalIndex)
##	dimnames(viterbiResults) <- list(featureNames(object), sampleNames(object))
##	if(returnSegments){
##		viterbiResults <- breaks(x=viterbiResults,
##					 states=states,
##					 position=position(object),
##					 chromosome=chromosome(object),
##					 verbose=verbose)
##	}
##	return(viterbiResults)
##}

setMethod("update", "environment", function(object, ...){
	if(length(ls(object)) == 0) stop("nothing to update")
	hmm(object=object[["locusset"]],
	    states=object[["states"]],
	    mu=object[["mu"]],
	    probs=object[["probs"]],
	    takeLog=object[["takeLog"]],
	    initialP=object[["initialP"]],
	    returnSegments=object[["returnSegments"]],
	    TAUP=object[["TAUP"]],
	    verbose=object[["verbose"]],
	    ice=object[["ice"]],
	    envir=object)
})

findFatherMother <- function(offspringId, object){
	stopifnot(!missing(offspringId))
	family.id <- pData(object)[sampleNames(object) == offspringId, "familyId"]
	father.id <- pData(object)[sampleNames(object) == offspringId, "fatherId"]
	mother.id <- pData(object)[sampleNames(object) == offspringId, "motherId"]
	father.name <- sampleNames(object)[object$familyId == family.id & object$individualId == father.id]
	mother.name <- sampleNames(object)[object$familyId == family.id & object$individualId == mother.id]
	if(length(father.name) > 1 | length(mother.name) > 1){
		stop("More than 1 father and/or more than 1 mother.  Check annotation in phenoData")
	}
	if(length(father.name) < 1 ){
		father.name <- NA
	}
	if(length(mother.name) < 1){
		mother.name <- NA
	}
	fmo.trio <- c(father.name, mother.name, offspringId)
	names(fmo.trio) <- c("father", "mother", "offspring")
	return(fmo.trio)
}
