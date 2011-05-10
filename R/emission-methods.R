##setMethod("calculateEmission", "oligoSnpSet", function(object,
calculateEmission.copynumber <- function(object, hmmOptions){
	cnStates <- hmmOptions[["copynumberStates"]]
	verbose <- hmmOptions[["verbose"]]
	states <- hmmOptions[["states"]]
	fn <- featureNames(object)
	S <- length(states)
	CN <- assayData(object)[["copyNumber"]]
	if(any(rowSums(is.na(CN)) == ncol(CN))){
		stop("Some rows have all missing values. Exclude these before continuing.")
	}
	if(any(colSums(is.na(CN)) == nrow(CN))){
		stop("Some samples have all missing values. Exclude these samples before continuing.")
	}
	sds <- 1/cnConfidence(object)
	tmp <- rowSums(!is.finite(sds))
	if(any(sds == 0, na.rm=TRUE)){
		warning("some sds were zero.  Replacing with typical value")
		sds[sds == 0] <- median(sds, na.rm=TRUE)
	}
	emission.cn <- array(NA, dim=c(nrow(object), ncol(object), S))
	if(!is.matrix(cnStates))
		cnStates <- matrix(cnStates, nrow(object), length(cnStates), byrow=TRUE)
	for(j in 1:ncol(object)){
		cn <- matrix(CN[, j], nrow(object), ncol(cnStates))
		sd <- matrix(sds[, j], nrow(object), ncol(cnStates))
		k <- which(!is.na(as.numeric(cn)))
		##emission.cn <- rep(NA, length(as.vector(cnStates)))
		tmp <- rep(NA, length(as.numeric(cnStates)))
		tmp[k] <- dnorm(x=as.numeric(cn)[k],
				mean=as.numeric(cnStates)[k],
				sd=as.numeric(sd)[k])
		emission.cn[, j, ] <- tmp
	}
	dimnames(emission.cn) <- list(featureNames(object),
				      sampleNames(object),
				      states)
	return(log(emission.cn))
}

calculateEmission.copynumber2 <- function(object, hmmOptions, prOutlier=0.01){
	cnStates <- hmmOptions[["copynumberStates"]]
	verbose <- hmmOptions[["verbose"]]
	states <- hmmOptions[["states"]]
	is.log <- hmmOptions[["is.log"]]
	fn <- featureNames(object)
	S <- length(states)
	CN <- assayData(object)[["copyNumber"]]
	if(any(rowSums(is.na(CN)) == ncol(CN))){
		stop("Some rows have all missing values. Exclude these before continuing.")
	}
	if(any(colSums(is.na(CN)) == nrow(CN))){
		stop("Some samples have all missing values. Exclude these samples before continuing.")
	}
	sds <- 1/cnConfidence(object)
	tmp <- rowSums(!is.finite(sds))
	if(any(sds == 0, na.rm=TRUE)){
		warning("some sds were zero.  Replacing with typical value")
		sds[sds == 0] <- median(sds, na.rm=TRUE)
	}
	emission.cn <- array(NA, dim=c(nrow(object), ncol(object), S))
	if(!is.matrix(cnStates))
		cnStates <- matrix(cnStates, nrow(object), length(cnStates), byrow=TRUE)
	if(is.log){
		MIN.CN <- -10
		MAX.CN <- 2.5
	} else {
		MIN.CN <- 0
		MAX.CN <- 10
	}
	if(any(CN < MIN.CN, na.rm=TRUE)) CN[CN < MIN.CN] <- MIN.CN
	if(any(CN > MAX.CN, na.rm=TRUE)) CN[CN > MAX.CN] <- MAX.CN
	for(j in 1:ncol(object)){
		cn <- matrix(CN[, j], nrow(object), ncol(cnStates))
		sd <- matrix(sds[, j], nrow(object), ncol(cnStates))
		k <- which(!is.na(as.numeric(cn)))
		##emission.cn <- rep(NA, length(as.vector(cnStates)))
		old.tmp <- tmp <- rep(NA, length(as.numeric(cnStates)))
		cnvector <- as.numeric(cn)[k]
		tmp[k] <- (1-prOutlier) * dnorm(x=cnvector,
						mean=as.numeric(cnStates)[k],
						sd=as.numeric(sd)[k]) +
			   prOutlier * dunif(cnvector, MIN.CN, MAX.CN)
##		old.tmp[k] <- dnorm(x=cnvector,
##				    mean=as.numeric(cnStates)[k],
##				    sd=as.numeric(sd)[k])
		emission.cn[, j, ] <- tmp
	}
	dimnames(emission.cn) <- list(featureNames(object),
				      sampleNames(object),
				      states)
	return(log(emission.cn))
}

calculateEmission.CopyNumberSet <- function(object, hmmOptions){
	EMIT.THR <- hmmOptions[["EMIT.THR"]]
	states <- hmmOptions[["states"]]
	verbose <- hmmOptions[["verbose"]]
	normalIndex <- hmmOptions[["normalIndex"]]
	if(all(is.na(cnConfidence(object)))){
		message("cnConfidence missing.  Using MAD")
		sds <- robustSds(copyNumber(object))
		cnConfidence(object) <- 1/sds
	}
	log.emission <- calculateEmission.copynumber2(object,
						      hmmOptions)
	if(any(is.na(log.emission))){
		if(verbose) message("Converting missing values in the emission matrix to 0")
		log.emission[is.na(log.emission)] <- 0
	}
##	if(any(log.emission < EMIT.THR)){
##		if(verbose) message("Minimum value in log emission probabilities is ", EMIT.THR, ".  See EMIT.THR in hmmOptions.")
##		log.emission[log.emission < EMIT.THR] <- EMIT.THR
##	}
	hmmOptions[["log.emission"]] <- log.emission
	hmmOptions
}


calculateEmission.SnpSet <- function(object, hmmOptions){
	ICE <- hmmOptions[["ICE"]]
	EMIT.THR <- hmmOptions[["EMIT.THR"]]
	states <- hmmOptions[["states"]]
	verbose <- hmmOptions[["verbose"]]
	normalIndex <- hmmOptions[["normalIndex"]]
	if(!ICE){
		log.gt.emission <- calculateEmission.genotype(object, hmmOptions)
	} else {
		##assumed order
		## ROH, normal
		log.gt.emission <- array(NA, dim=c(nrow(object), ncol(object), length(states)),
					 dimnames=list(featureNames(object),
					 sampleNames(object),
					 states))
		tmp <- genotypeEmissionCrlmm(object, hmmOptions)
		rohStates <- which(hmmOptions[["rohStates"]])
		notRohState <- which(!hmmOptions[["rohStates"]])
		for(j in rohStates){
			log.gt.emission[, , j] <- tmp[, , "ROH"]
		}
		for(j in notRohState){
			log.gt.emission[, , j] <- tmp[, , "normal"]
		}
	}
	log.emission <- log.gt.emission
	rm(log.gt.emission); gc()
	if(any(is.na(log.emission))){
		if(verbose) message("Converting missing values in the emission matrix to 0")
		log.emission[is.na(log.emission)] <- 0
	}
	if(any(log.emission < EMIT.THR)){
		if(verbose) message("Minimum value in log emission probabilities is ", EMIT.THR, ".  See EMIT.THR in hmmOptions.")
		log.emission[log.emission < EMIT.THR] <- EMIT.THR
	}
	hmmOptions[["log.emission"]] <- log.emission
	hmmOptions
}


calculateEmission.oligoSnpSet <- function(object, hmmOptions){
	ICE <- hmmOptions[["ICE"]]
	EMIT.THR <- hmmOptions[["EMIT.THR"]]
	states <- hmmOptions[["states"]]
	verbose <- hmmOptions[["verbose"]]
	normalIndex <- hmmOptions[["normalIndex"]]
	if(all(is.na(cnConfidence(object)))){
		message("cnConfidence missing.  Using MAD")
		sds <- robustSds(copyNumber(object))
		cnConfidence(object) <- 1/sds
	} ##else {
	##sds <- 1/cnConfidence(object)
##	}
	log.cn.emission <- calculateEmission.copynumber2(object,
							 hmmOptions)
	if(!ICE){
		log.gt.emission <- calculateEmission.genotype(object, hmmOptions)
	} else {
		##assumed order
		## ROH, normal
		log.gt.emission <- array(NA, dim=c(nrow(object), ncol(object), length(states)),
					 dimnames=list(featureNames(object),
					 sampleNames(object),
					 states))
		tmp <- genotypeEmissionCrlmm(object, hmmOptions)
		rohStates <- which(hmmOptions[["rohStates"]])
		notRohState <- which(!hmmOptions[["rohStates"]])
		for(j in rohStates){
			log.gt.emission[, , j] <- tmp[, , "ROH"]
		}
		for(j in notRohState){
			log.gt.emission[, , j] <- tmp[, , "normal"]
		}
	}
	log.emission <- log.gt.emission+log.cn.emission
	if(any(is.na(log.emission))){
		if(verbose) message("Converting missing values in the emission matrix to 0")
		log.emission[is.na(log.emission)] <- 0
	}
##	if(any(log.emission < EMIT.THR)){
##		if(verbose) message("Minimum value in log emission probabilities is ", EMIT.THR, ".  See EMIT.THR in hmmOptions.")
##		log.emission[log.emission < EMIT.THR] <- EMIT.THR
##	}
	hmmOptions[["log.gt.emission"]] <- log.gt.emission
	hmmOptions[["log.cn.emission"]] <- log.cn.emission
	hmmOptions[["log.emission"]] <- log.emission
	hmmOptions
}

setMethod("calculateEmission", "CopyNumberSet",
	  function(object, hmmOptions){
		  calculateEmission.CopyNumberSet(object, hmmOptions)
	  })


setMethod("calculateEmission", "oligoSnpSet",
	  function(object, hmmOptions){
		  calculateEmission.oligoSnpSet(object, hmmOptions)
	  })

setMethod("calculateEmission", "SnpSet",
	  function(object, hmmOptions){
		  calculateEmission.SnpSet(object, hmmOptions)
	  })

##setMethod("calculateEmission", "SnpCopyNumberSet", function(object, mu, states, envir, ...){
##	if(all(is.na(cnConfidence(object)))){
##		sds <- robustSds(copyNumber(object))
##	} else {
##		sds <- 1/cnConfidence(object)
##	}
##	emission.cn <- copynumberEmission(copynumber=copyNumber(object),
##					  mu=mu,
##					  sds=sds,
##					  states=states,
##					  takeLog=FALSE,
##					  verbose=verbose)
##	envir[["emission.cn"]] <- emission.cn
##})

##setMethod("calculateEmission", "SnpCallSet", function(object, probs, probMissing, states, envir, ice=FALSE){
##	if(!ice){
##		emission.gt <- genotypeEmission(genotypes=calls(object),
##						probHomCall=probs,
##						probMissing=probMissing,
##						states=states)
##	} else{
##		emission.gt <- genotypeEmissionCrlmm(genotypes=calls(object),
##						     conf=callsConfidence(object),
##						     pHetCalledHom=probs["pHetCalledHom"],
##						     pHetCalledHet=probs["pHetCalledHet"],
##						     pHomInNormal=probs["pHomInNormal"],
##						     pHomInRoh=probs["pHomInRoh"])
##	}
##	envir[["emission.gt"]] <- emission.gt
##})

calculateEmission.genotype <- function(object, hmmOptions){
	states <- hmmOptions[["states"]]
	p <- hmmOptions[["prGenotypeHomozygous"]]
	prGenotypeMissing <- hmmOptions[["prGenotypeMissing"]]
	verbose <- hmmOptions[["verbose"]]
	stopifnot(length(p) == length(states))
	if(!is.numeric(calls(object))) stop("genotypes must be integers (1=AA, 2=AB, 3=BB) or NA (missing)")
	GT <- calls(object)
	emission <- array(GT, dim=c(nrow(GT), ncol(GT), length(states)), dimnames=list(featureNames(object), sampleNames(object), states))
	missingGT <- any(is.na(GT))
	if(missingGT){
		if(verbose) message("Some genotypes are NAs.  The default assumes that prGenotypeMissing is the same for each state -- see hmmOptions")
	}
	for(s in seq(along=states)){
		tmp <- GT
		tmp[tmp == 1 | tmp == 3] <- p[s]
		tmp[tmp == 2] <- 1-p[s]
		if(missingGT){
			tmp[is.na(tmp)] <- prGenotypeMissing[s]
		}
		emission[, , s] <- tmp
	}
	logemit <- log(emission)
	return(logemit)
}

##genotypeEmission <- function(genotypes, conf, states, probHomCall,
##			     probMissing, verbose=TRUE){
##	##function to replace .getCallEmission
##	##message("Use genotypeEmissionCrlmm if using crlmm-processed data on affy250k or affy 6.0")
##	if(!is.numeric(genotypes)) stop("genotypes must be integers (1=AA, 2=AB, 3=BB, 4=missing")
##	emission <- array(genotypes, dim=c(nrow(genotypes), ncol(genotypes), length(states)))
##	for(s in seq(along=states)){
##		tmp <- genotypes
##		tmp[tmp == 1 | tmp == 3] <- probHomCall[s]
##		tmp[tmp == 2] <- 1-probHomCall[s]
##		if(!missing(probMissing)) tmp[tmp == 4 | is.na(tmp)] <- probMissing[s]
##		emission[, , s] <- tmp
##	}
##	dimnames(emission) <- list(rownames(genotypes),
##				   colnames(genotypes),
##				   states)
##	logemit <- log(emission)
##	return(logemit)
##}


icePlatforms <- function(){
	c("pd.genomewidesnp.6",
	  "genomewidesnp6",
	  "pd.mapping250k.nsp",
	  "pd.mapping250k.sty",
	  "pd.mapping250k.nsp, pd.mapping250k.sty")
}

##genotypeEmissionCrlmm <- function(genotypes, conf,
##				  pHetCalledHom=0.001,
##				  pHetCalledHet=0.995,
##				  pHomInNormal=0.8,
##				  pHomInRoh=0.999,  ##Pr(AA or BB | region of homozygosity)
##				  annotation){
genotypeEmissionCrlmm <- function(object, hmmOptions){
	if(!annotation(object) %in% icePlatforms()){
		message("ICE is TRUE, but hapmap crlmm confidence scores for ", annotation(object), " are not available. Using crlmm confidence scores from HapMap samples assayed on the Affy 6.0 platform.")
		annotation <- "genomewidesnp6"
	} else {
		if(annotation(object) == "pd.genomewidesnp.6"){
			annotation <- "genomewidesnp6"
		} else annotation <- annotation(object)
	}
	loader(paste(annotation, "Conf.rda", sep=""), .vanillaIcePkgEnv, "VanillaICE")
	hapmapP <- getVarInEnv("reference")
	pHetCalledHom <- hmmOptions[["pHetCalledHom"]]
	pHetCalledHet <- hmmOptions[["pHetCalledHet"]]
	pHomInNormal <- hmmOptions[["pHomInNormal"]]
	pHomInRoh <- hmmOptions[["pHomInRoh"]]
	if(length(annotation(object)) < 1) stop("must specify annotation")
	GT <- as.integer(calls(object))
	GTconf <- confs(object)
	##data(list=paste(annotation, "Conf", sep=""), package="VanillaICE", envir=environment())
	if(length(pHomInNormal) == nrow(GTconf)){  ##convert to vector
		pHomInNormal <- as.numeric(matrix(pHomInNormal, nrow(GTconf), ncol(GTconf), byrow=FALSE))
	} else pHomInNormal <-  rep(pHomInNormal, length(GT))
	hapmapP[, 2] <- 1-exp(-hapmapP[, 2]/1000)
	##p = 1-exp(-X/1000)
	##1000*log(1-p)=X
	##confidence <- 1-exp(-GTconf/1000)
	i11 <- hapmapP[, 1] == 3  ##called homozygous truth homozygous
	i12 <- hapmapP[, 1] == 4  ##called homozygous truth heterozygous
	i21 <- hapmapP[, 1] == 1  ##called het truth hom
	i22 <- hapmapP[, 1] == 2  ##called het truth het
	f11 <- density(hapmapP[i11, 2], from=0, to=1, n=1e3)
	f12 <- density(hapmapP[i12, 2], from=0, to=1, n=1e3)
	f21 <- density(hapmapP[i21, 2], from=0, to=1, n=1e3)
	f22 <- density(hapmapP[i22, 2], from=0, to=1, n=1e3)
	##-------------------------------------------------------------------------
	##distribution of observed call probabilities when the call is homozygous
	##-------------------------------------------------------------------------
	##-------------------------------------------------------------------------
	##P(phat | LOH, gte)
	##-------------------------------------------------------------------------
	##GT <- as.integer(genotypes)
	##confidence <- as.numeric(confidence)
	confidence <- as.numeric(GTconf)
	pTruthIsNormal <- pTruthIsRoh <- rep(NA, length(GT))
	confidence[confidence==0] <- 0.01 ##Otherwise, NA's result
	hom <- which(GT == 1 | GT == 3)
	observedPcalledHom <- cut(confidence[hom], breaks=f11$x, labels=FALSE)
	pTruthIsRoh[hom] <- f11$y[observedPcalledHom]
	het <- which(GT == 2)
	observedPcalledHet <- cut(confidence[het], breaks=f11$x, labels=FALSE)
	pTruthIsRoh[het] <- f21$y[observedPcalledHet]
	##-------------------------------------------------------------------------
	##Calculate P(phat | Normal, HOM)
	##-------------------------------------------------------------------------
	chet1 <- f22$y[cut(confidence[het], breaks=f22$x, labels=FALSE)]
	chet2 <- f21$y[cut(confidence[het], breaks=f21$x, labels=FALSE)]
	##term5[1]=P(true genotype is HET | genotype call is AB, state is normal)
	pTruthIsNormal[het] <- chet1*pHetCalledHet + chet2*(1-pHetCalledHet)
	##chom1=called homozygous truth heterozygous
	chom1 <- f12$y[cut(confidence[hom], breaks=f12$x, labels=FALSE)]
	##chom2=called homozygous truth homozygous
	chom2 <- f11$y[cut(confidence[hom], breaks=f11$x, labels=FALSE)]
	##chom4 <- 0.9999    ##P(HOM|CHOM)
	##probability that the true state is HOM when genotype call is homozygous
	##pHetCalledHom = P(true genotype is HET | calls is AA or BB, state is normal)
	pTruthIsNormal[hom] <- chom1*pHetCalledHom + chom2*(1-pHetCalledHom)
	fNormal <- fLoh <- rep(NA, length(GT))
	fNormal[hom] <- pHomInNormal[hom] * pTruthIsNormal[hom]
	fNormal[het] <- (1-pHomInNormal[het]) * pTruthIsNormal[het]
	fLoh[hom] <- pHomInRoh * pTruthIsRoh[hom]
	fLoh[het] <- (1-pHomInRoh) * pTruthIsRoh[het]
	f <- array(NA, dim=c(nrow(object), ncol(object), 2), dimnames=list(featureNames(object),
							     sampleNames(object),
							     c("normal", "ROH")))
	f[, , "normal"] <- matrix(fNormal, nrow(object), ncol(object))
	f[, , "ROH"] <- matrix(fLoh, nrow(object), ncol(object))
	f[f  == 0] <- min(f[f > 0], na.rm=TRUE)
	f <- log(f)
	return(f)
}

hmm.SnpSuperSet <- function(object, hmmOptions){
	if(ncol(object) < 3){
		return("No complete trios")
	}
	stopifnot(all(c("familyId", "fatherId", "motherId", "individualId") %in% varLabels(object)))
	TAUP <- hmmOptions[["TAUP"]]
	states <- hmmOptions[["states"]]
	log.initial <- hmmOptions[["log.initial"]]
	verbose <- hmmOptions[["verbose"]]
	normal2altered <- hmmOptions[["normal2altered"]]
	altered2normal <- hmmOptions[["altered2normal"]]
	normalIndex <- hmmOptions[["normalIndex"]]
	trios <- hmmOptions[["trios"]]
	##For each offspring, find the father and mother in the same family
	rD <- vector("list", nrow(trios))
	for(i in 1:nrow(trios)){
		##if(verbose) cat("Family ", unique(familyId)[i], ", ")
		if(verbose) cat("Offspring ID ", trios[i, "offspring"], "\n")
		trioSet <- object[, match(trios[i, ], sampleNames(object))]
		## Remove the noinformative snps here.
		isBPI <- isBiparental.SnpSuperSet(trioSet, allowHetParent=hmmOptions[["allowHetParent"]])
		isInformative <- !is.na(isBPI)
		if(all(!isInformative)){
			fit[, i] <- 1
			next()
		}
		trioSet <- trioSet[isInformative, ]
		isBPI <- isBPI[isInformative]
		log.emission <- computeBpiEmission.SnpSuperSet(trioSet, hmmOptions, isBPI=isBPI)
		##.GlobalEnv[["emission"]] <- emission
		if(is.null(log.emission)) stop("not a father, mother, offspring trio")
		log.e <- array(log.emission, dim=c(nrow(trioSet), 1, 2), dimnames=list(featureNames(trioSet), trios[i, "offspring"], hmmOptions[["states"]]))
		hmmOptions[["log.emission"]] <- log.e
		index <- match(featureNames(trioSet), featureNames(object))
		rD[[i]] <- viterbi(trioSet, hmmOptions)
	}
	L <- sapply(rD, nrow)
	if(any(L == 1) & any(L > 1)){
		rD <- c(do.call(c, rD[L == 1]), do.call(c, rD[L > 1]))
	} else {
		rD <- do.call(c, rD)
	}
	return(rD)
}
##		vitResults <- viterbi(initialStateProbs=log(initialP),
##				      emission=log.e,
##				      tau=tau[, "transitionPr"],
##				      arm=tau[, "arm"],
##				      normalIndex=normalIndex,
##				      verbose=verbose,
##				      normal2altered=normal2altered,
##				      altered2normal=altered2normal,
##				      returnLikelihood=TRUE)
		##vitSequence is a vector -- one trio at a time
##		vitSequence <- vitResults[["stateSequence"]]
##		if(length(table(tau[, "arm"])) > 1){
##			##insert an extra index to force a break between chromosome arms
##			tmp <- rep(NA, length(vitSequence)+1)
##			end.parm <- end(Rle(tau[, "arm"]))[1]
##			tmp[1:end.parm] <- vitSequence[1:end.parm]
##			tmp[end.parm+1] <- 999
##			tmp[(end.parm+2):length(tmp)] <- vitSequence[((end.parm)+1):length(vitSequence)]
##			vitSequence <- tmp
##		}
##		llr <- vitResults[["logLikelihoodRatio"]][[1]]
##		rl <- Rle(vitSequence)
##		start.index <- start(rl)[runValue(rl) != 999]
##		end.index <- end(rl)[runValue(rl) != 999]
##		##this is tricky since we've added an index to force a segment for each arm.
##		armBreak <- which(vitSequence==999)
##		if(length(armBreak) > 0){
##			start.index[start.index > armBreak] <- start.index[start.index > armBreak] - 1
##			end.index[end.index > armBreak] <- end.index[end.index > armBreak] - 1
##		}
##		start <- position(trioSet)[start.index]
##		end <- position(trioSet)[end.index]
##		##numMarkers <- unlist(numMarkers)
##		numMarkers <- width(rl)[runValue(rl) != 999]
##		states <- hmmOptions[["states"]][vitSequence[start.index]]
##		##states <- (hmmOptions[["states"]])[vitSequence[start.index]]
##		ir <- IRanges(start=start, end=end)
##		##For each segment, calculate number biparental, number not biparental
##		nBpi <- nNotBpi <- rep(NA, length(ir))
##		for(j in 1:length(start)){
##			region <- (start(rl)[j]):(end(rl)[j])
##			region <- (start.index[j]):(end.index[j])
##			nNonInformative <- sum(is.na(isBPI[region]))
##			nInformative <- sum(!is.na(isBPI[region]))
##			nNotBpi[j] <- sum(isBPI[region] == FALSE, na.rm=TRUE)
##			nBpi[j] <- sum(isBPI[region] == TRUE, na.rm=TRUE)
##		}
##		rD[[i]] <- RangedData(ir,
##				      space=rep(paste("chr", unique(chromosome(trioSet)), sep=""), length(ir)),
##				      offspringId=rep(trios[i, "offspring"], length(ir)),
##				      state=states,
##				      numMarkers=numMarkers,
##				      nNotBpi=nNotBpi,
##				      nBpi=nBpi,
##				      LLR=llr)
##	}
	##to avoid a .Primivite error with do.call(c, rD)
##	tmp <- do.call(c, rD[sapply(rD, nrow) == 1])
##	tmp2 <- do.call(c, rD[sapply(rD, nrow) > 1])
##	rD <- c(tmp, tmp2)
##	return(rD)
#}

computeBpiEmission.SnpSuperSet <- function(object, hmmOptions, isBPI){
	states <- hmmOptions[["states"]]
	prGtError <- hmmOptions[["prGtError"]]
	ICE <- hmmOptions[["ICE"]]
	emission <- matrix(NA, nrow(object), ncol=2)
	colnames(emission) <- states
	if(ICE){
		pCrlmm <- confs(object)  ## crlmm confidence score
		## take the minimum confidence score in the trio
		pCrlmm <- apply(pCrlmm, 1, min, na.rm=TRUE)
		## set emission probability to min(crlmmConfidence, 0.999)
		I <- as.integer(pCrlmm < (1 - prGtError[["BPI"]]))
		pCrlmm <- pCrlmm*I + (1 - prGtError[["BPI"]])*(1-I)
		emission[,  "BPI"] <- pCrlmm
		##Pr(mendelian inconsistency | BPI) = 0.001
		emission[, "BPI"] <- 1 - pCrlmm
	} else { ##ignore confidence scores
		##Pr(call is consistent with biparental inheritance | BPI) = 0.999
		emission[isBPI==TRUE,  "BPI"] <-  1-prGtError["BPI"]
		##Pr(mendelian inconsistency | BPI) = 0.001
		emission[isBPI==FALSE, "BPI"] <- prGtError["BPI"] ##Mendelian inconsistancy
	}
	##Pr(call is consistent with biparental inheritance | not BPI) = 0.01
	emission[isBPI==TRUE,  "notBPI"] <- prGtError["notBPI"]   ## biparental inheritance, but true state is not Biparental
	##Pr(mendelian inconsistency | not BPI) = 0.99
	emission[isBPI==FALSE, "notBPI"] <- 1-prGtError["notBPI"] ## Mendelian inconsistancy
	log.emission <- log(emission)
	return(log.emission)
}




isBiparental.SnpSuperSet <- function(object, allowHetParent=FALSE){
	##if(length(object$familyMember) < 3) stop("object$familyMember not the right length")
	father <- 1
	mother <- 2
	offspring <- 3
	F <- calls(object[, father])
	M <- calls(object[, mother])
	O <- calls(object[, offspring])
	object <- cbind(F, M, O)
	colnames(object) <- c("father", "mother", "offspring")
	biparental <- isBiparental.matrix(object, allowHetParent=allowHetParent)
	return(biparental)
}

isBiparental.matrix <- function(object, allowHetParent=TRUE){
	F <- object[, 1]
	M <- object[, 2]
	O <- object[, 3]
	##M/F AA, F/M BB, O AB
	##isHet <- offspringHeterozygous(object)  ##offspring is heterozygous
	biparental <- rep(NA, nrow(object))
	biparental[F==1 & M == 3 & O == 2] <- TRUE#Pr(O | biparental)=1-0.001, Pr(O | not biparental) = 0.01
	biparental[F==3 & M == 1 & O == 2] <- TRUE#Pr(O | biparental)=1-0.001, Pr(O | not biparental) = 0.01
	##M/F AA, F/M BB, O AA or BB
	biparental[F==1 & M == 3 & (O == 1 | O == 3)] <- FALSE#Pr(O | biparental)=0.001, Pr(O | not biparental) = 1-0.01
	biparental[F==3 & M == 1 & (O == 1 | O == 3)] <- FALSE#Pr(O | biparental)=0.001, Pr(O | not biparental) = 1-0.01
	## M/F AA, F/M BB, O AB
	if(allowHetParent) biparental[F == 1 & M == 2 & O == 2] <- TRUE#Pr(O | biparental)=1-0.001, Pr(O | not biparental) = 0.01
	if(allowHetParent) biparental[F == 2 & M == 1 & O == 2] <- TRUE#Pr(O | biparental)=1-0.001, Pr(O | not biparental) = 0.01
	## F AB, M AA, O BB is not biparental
	## F AA, M AB, O BB is not biparental
	biparental[F == 2 & M == 1 & O == 3] <- FALSE#Pr(O | biparental)=1-0.001, Pr(O | not biparental) = 0.01
	biparental[F == 1 & M == 2 & O == 3] <- FALSE#Pr(O | biparental)=0.001, Pr(O | not biparental) = 1-0.01
	## M AA, F AB, O AB
	if(allowHetParent) biparental[F == 2 & M == 3 & O == 2] <- TRUE#Pr(O | biparental)=1-0.001, Pr(O | not biparental) = 0.01
	if(allowHetParent) biparental[F == 3 & M == 2 & O == 2] <- TRUE#Pr(O | biparental)=1-0.001, Pr(O | not biparental) = 0.01
	## F=AB, M=BB, O=AA is NOT biparental
	biparental[F == 2 & M == 3 & O == 1] <- FALSE#Pr(O | biparental)=0.001, Pr(O | not biparental) = 1-0.01
	biparental[F == 3 & M == 2 & O == 1] <- FALSE#Pr(O | biparental)=0.001, Pr(O | not biparental) = 1-0.01
	return(biparental)
}

##computeHmm.CNSet <- function(object, cnOptions){
##	hmmOptions <- cnOptions[["hmmOpts"]]
##	object <- object[order(chromosome(object), position(object)), ]
##	##emission <- hmmOptions[["emission"]]
##	chrom <- unique(chromosome(object))
####	tPr <- transitionProbability(chromosome=chromosome(object),
####				     position=position(object),
####				     TAUP=hmmOptions[["TAUP"]])
####	emissionPr(object) <- computeEmission(object, hmmOptions)
##	rangedData(object) <- viterbi.CNSet(object,
##					    hmmOptions=hmmOptions,
##					    transitionPr=tPr[, "transitionPr"],
##					    chromosomeArm=tPr[, "arm"])
##	return(object)
##}

##viterbi.CNSet <- function(object, hmmOptions, transitionPr, chromosomeArm){
##	state.sequence <- viterbi(object,
##				  hmmOptions)
####				  emission=emissionPr(object),
######				  tau=transitionPr,
####				  initialStateProbs=hmmOptions[["log.initial"]],
####				  arm=chromosomeArm,
####				  normalIndex=hmmOptions[["normalIndex"]],
####				  normal2altered=hmmOptions[["normal2altered"]],
####				  altered2normal=hmmOptions[["altered2normal"]],
####				  altered2altered=hmmOptions[["altered2altered"]])
##	state.sequence <- data.frame(state.sequence)
##	rleList <- RleList(state.sequence)
##	rd <- RangedData(rleList)
####	rdList <- vector("list", length(rle.object))
####	for(i in seq(along=rdList)){
####		rdList[[i]] <- RangedData(rle.object[[i]],
####					  space=paste("chr", transitionPr[, "chromosome"], sep=""),
####					  state=runValue(rle.object[[i]]),
####					  sample=sampleNames(object)[i],
####					  nprobes=runLength(rle.object[[i]]))
####	}
####	rangedData <- do.call("c", rdList)
##	return(rd)
##}

calculateEmission.CNSet <- function(object, hmmOptions){
	EMIT.THR <- hmmOptions[["EMIT.THR"]]
	cnStates <- hmmOptions[["copynumberStates"]]
	object <- object[order(chromosome(object), position(object)), ]
	if(any(diff(position(object)) < 0)) stop("must be ordered by chromosome and physical position")
	emissionProbs <- array(NA, dim=c(nrow(object), ncol(object), length(hmmOptions[["copynumberStates"]])))
	dimnames(emissionProbs) <- list(featureNames(object),
					sampleNames(object),
					paste("copy.number_", hmmOptions[["copynumberStates"]], sep=""))
	batch <- object$batch
	for(i in seq(along=unique(batch))){
		emissionProbs[, batch == unique(batch)[i], ] <- getEmission(object[, batch==unique(batch)[i]], hmmOptions)
	}
	if(EMIT.THR > -Inf){  ## truncate emission probabilities for outliers
		emissionProbs[emissionProbs < EMIT.THR] <- EMIT.THR
	}
	hmmOptions[["log.emission"]] <- emissionProbs
	hmmOptions
}

getEmission <- function(object, hmmOptions){
	emissionProbs <- array(NA, dim=c(nrow(object),
				   ncol(object), length(hmmOptions[["copynumberStates"]])))
	emit.snps <- getEmission.snps(object[isSnp(object), ], hmmOptions)
	if(any(!isSnp(object))){
		emit.nps <- getEmission.nps(object[!isSnp(object), ], hmmOptions)
		emissionProbs[!isSnp(object), , ] <- emit.nps
	}
	emissionProbs[isSnp(object), , ] <- emit.snps
	emissionProbs
}

getEmission.nps <- function(object, hmmOptions){
	##****************************************************
	##	                                             *
	##  Emission probabilities for nonpolymorphic probes *
	##	                                             *
	##****************************************************
	batch <- unique(object$batch)
	scaleSds <- hmmOptions[["scaleSds"]]
	cnStates <- hmmOptions[["copynumberStates"]]
	verbose <- hmmOptions[["verbose"]]
	if(verbose) message("Computing emission probabilities for nonpolymorphic loci.")
	if(scaleSds){
		##a <- log2(CA(object))
		##sds.a <- apply(a, 2, mad, na.rm=TRUE)
		##sds.a <- sds.a/median(sds.a)
		sds.a <- robustSds(log2(CA(object)))
		##sds.a[sds.a < 1] <- 1
		sds.a <- matrix(sds.a, nrow(object), ncol(object), byrow=TRUE)
	} else sds.a <- matrix(0, nrow(object), ncol(object))
	emissionProbs <- array(NA, dim=c(nrow(object),
				   ncol(object), length(cnStates)))
	nuA <- getParam(object, "nuA", batch)
	phiA <- getParam(object, "phiA", batch)
	sig2A <- getParam(object, "sig2A", batch)
	if(any(is.na(sig2A))){
		sig2A[is.na(sig2A)] <- median(sig2A, na.rm=TRUE)
	}
	##tau2A <- getParam(object, "tau2A", batch)
	##Assume that on the log-scale, that the background variance is the same...
	##tau2A <- sig2A
	a <- as.numeric(log2(A(object)))
##	if(any(cnStates > 2)){
##		cnStates[cnStates > 2] <- cnStates[cnStates > 2] * 0.85
##	}
	for(k in seq(along=cnStates)){
		CT <- cnStates[k]
		mus.matrix=matrix(log2(nuA + CT*phiA), nrow(object), ncol(object))
		mus <- as.numeric(matrix(log2(nuA + CT*phiA), nrow(object), ncol(object)))
		sds.matrix <- matrix(sqrt(sig2A), nrow(object), ncol(object))
		##sds.matrix <- sds.matrix + sds.a
		sds.matrix <- sds.matrix*sds.a
		sds <- as.numeric(sds.matrix)
		tmp <- matrix(dnorm(a, mean=mus, sd=sds), nrow(object), ncol(object))
		emissionProbs[, , k] <- log(tmp)
	}
	emissionProbs
}


getEmission.snps <- function(object, hmmOptions){
	batch <- unique(object$batch)
	if(length(batch) > 1) stop("batch variable not unique")
	scaleSds <- hmmOptions[["scaleSds"]]
	cnStates <- hmmOptions[["copynumberStates"]]
	verbose <- hmmOptions[["verbose"]]
	if(verbose) message("Computing emission probabilities for polymorphic loci.")
	if(scaleSds){
		##a <- log2(CA(object) + CB(object))
		##sds.a <- apply(a, 2, mad, na.rm=TRUE)
		##sds.a <- sds.a/median(sds.a)
		sds.a <- robustSds(log2(CA(object) + CB(object)))
		##sds.a[sds.a < 1] <- 1
		sds.a <- matrix(sds.a, nrow(object), ncol(object), byrow=TRUE)
	} else sds.a <- sds.b <- matrix(0, nrow(object), ncol(object))
	emissionProbs <- array(NA, dim=c(nrow(object),
				   ncol(object), length(cnStates)))
	corr <- getParam(object, "corr", batch)
	corrA.BB <- getParam(object, "corrA.BB", batch)
	corrB.AA <- getParam(object, "corrB.AA", batch)
	nuA <- getParam(object, "nuA", batch)
	nuB <- getParam(object, "nuB", batch)
	phiA <- getParam(object, "phiA", batch)
	phiB <- getParam(object, "phiB", batch)
	sig2A <- getParam(object, "sig2A", batch)
	sig2B <- getParam(object, "sig2B", batch)
	tau2A <- getParam(object, "tau2A", batch)
	tau2B <- getParam(object, "tau2B", batch)
	a <- as.numeric(log2(A(object)))
	b <- as.numeric(log2(B(object)))
	for(k in seq(along=cnStates)){
		T <- cnStates[k]
		f.x.y <- matrix(0, sum(nrow(object)), ncol(object))
		for(copyA in 0:T){
			copyB <- T-copyA
			sigmaA <- sqrt(tau2A*(copyA==0) + sig2A*(copyA > 0))
			sigmaB <- sqrt(tau2B*(copyB==0) + sig2B*(copyB > 0))
			if(copyA == 0 & copyB > 0) r <- corrA.BB
			if(copyA > 0 & copyB == 0) r <- corrB.AA
			if(copyA > 0 & copyB > 0) r <- corr
			if(copyA == 0 & copyB == 0) r <- 0
			muA <- log2(nuA+copyA*phiA)
			muB <- log2(nuB+copyB*phiB)

			sigmaA <- matrix(sigmaA, nrow=length(sigmaA), ncol=ncol(object), byrow=FALSE)
			sigmaB <- matrix(sigmaB, nrow=length(sigmaB), ncol=ncol(object), byrow=FALSE)
			## scale the variances by a sample-specific estimate of the variances
			## var(I_A, ijp) = sigma_A_ip * sigma_A_jp
			##sigmaA <- sigmaA+sds.a
			##sigmaB <- sigmaB+sds.a
			sigmaA <- sigmaA*sds.a
			sigmaB <- sigmaB*sds.a
			meanA <- as.numeric(matrix(muA, nrow(object), ncol(object)))
			meanB <- as.numeric(matrix(muB, nrow(object), ncol(object)))
			rho <- as.numeric(matrix(r, nrow(object), ncol(object)))
			sd.A <- as.numeric(matrix(sigmaA, nrow(object), ncol(object)))
			sd.B <- as.numeric(matrix(sigmaB, nrow(object), ncol(object)))
			Q.x.y <- 1/(1-rho^2)*(((a - meanA)/sd.A)^2 + ((b - meanB)/sd.B)^2 - 2*rho*((a - meanA)*(b - meanB))/(sd.A*sd.B))
			## For CN states > 1, assume that any of the possible genotypes are equally likely a priori...just take the sum
			## For instance, for state copy number 2 there are three combinations: AA, AB, BB
			##   -- two of the three combinations should be near zero.
			## TODO: copy-neutral LOH would put near-zero mass on both copyA > 0, copyB > 0
			f.x.y <- f.x.y + matrix(1/(2*pi*sd.A*sd.B*sqrt(1-rho^2))*exp(-0.5*Q.x.y), nrow(object), ncol(object))
		}
		emissionProbs[, , k] <- log(f.x.y)
	}
	emissionProbs
}


setMethod("calculateEmission", "CNSet", function(object, hmmOptions){
	calculateEmission.CNSet(object, hmmOptions)
})

##setMethod("computeEmission", "character", function(object, hmmOptions){
##	filename <- object
##	chrom <- gsub(".rda", "", strsplit(filename, "_")[[1]][[2]])
##	if(hmmOptions[["verbose"]])
##		message("Compute emission probabilities for chromosome ", chrom)
##	if(file.exists(filename)){
##		load(filename)
##		cnSet <- get("cnSet")
##	} else {
##		stop("File ", filename, " does not exist.")
##	}
##	emission <- computeEmission(cnSet, hmmOptions)
##	message("Saving ", file.path(ldPath(), paste("emission_", chrom, ".rda", sep="")))
##	if(hmmOptions[["save.it"]]){
##		save(emission,
##		     file=file.path(ldPath(), paste("emission_", chrom, ".rda", sep="")))
##	}
##	return(emission)
##})

##setMethod("computeHmm", "CNSet", function(object, hmmOptions){
##	computeHmm.CNSet(object, hmmOptions)
##})
