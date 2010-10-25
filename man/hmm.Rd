\name{hmm}
\alias{hmm}
\alias{viterbi}

\title{Application of the Viterbi algorithm to copy number and/or
genotype data.}
\description{
  A wrapper for fitting the HMM.
}
\usage{
hmm(object, hmm.params, ...)
viterbi(object, hmm.params, verbose=TRUE, normal2altered=1, altered2normal=1,
altered2altered=1, TAUP=1e8)
}

\arguments{

\item{object}{one of the following classes derived from eSet:
	\code{SnpSet}, \code{oligoSnpSet}, \code{CopyNumberSet},
	\code{CNSet}}

\item{hmm.params}{List. See \code{hmm.setup}}

\item{verbose}{Logical.  Whether to display all messages and warnings.}

\item{normal2altered}{Numeric. Factor for scaling the probability of transitioning
	      from the normal state to an altered state}

\item{altered2normal}{Numeric. Factor for scaling the probability of transitioning
      from an altered state to a normal state.}

\item{altered2altered}{Numeric. Factor for scaling the probability of transitioning
      from an altered state to a different altered state}

\item{TAUP}{

  Numeric. Factor for scaling the distance weighted transition
  probability.  The transition probability is computed as

  \code{exp(-2 * d/TAUP)}, where d is the distance between two loci.

  The default value is 1e8, but users can adjust this number to
  control the smoothness of the resulting HMM.

}

  \item{\dots}{Additional arguments to \code{viterbi}.}

}

\details{

}

\value{

	An object of class \code{RangedData}.

}

\seealso{\code{\link{hmm.setup}}}

\references{

    RB Scharpf et al. (2008) Hidden Markov Models for the assessment of
  chromosomal alterations using high-throughput SNP arrays, Annals of
  Applied Statistics

  }

\examples{

data(locusLevelData)
oligoSet <- new("oligoSnpSet",
		copyNumber=log2(locusLevelData[["copynumber"]]/100),
		call=locusLevelData[["genotypes"]],
		callProbability=locusLevelData[["crlmmConfidence"]],
		annotation=locusLevelData[["platform"]])
oligoSet <- oligoSet[!is.na(chromosome(oligoSet)), ]
oligoSet <- oligoSet[order(chromosome(oligoSet), position(oligoSet)),]
hmmOpts <- hmm.setup(oligoSet,
		     copynumberStates=log2(c(1, 2, 2, 3)),
		     states=c("hem-del", "ROH", "normal", "amp"),
		     normalIndex=3,
		     log.initialP=rep(log(1/4), 4),
		     prGenotypeHomozygous=c(0.99, 0.99, 0.7, 0.7))
fit <- hmm(oligoSet, hmmOpts, TAUP=5e7)

## Useful accessors for RangedData
tmp <- ranges(fit)
start(fit)
end(fit)

##Log likelihood ratio comparing likelihood of predicted state to the 'normal' state
## for each segment
fit$LLR
## the number of SNPs / nonpolymorphic loci in each segment
fit$numMarkers

## Expand to a locus x sample matrix


}

\author{R. Scharpf}

\keyword{models}
\keyword{manip}
\keyword{ts}

