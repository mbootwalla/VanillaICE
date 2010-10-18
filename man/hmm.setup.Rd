\name{hmm.setup}
\alias{hmm.setup}
\title{
	Wrapper for computing emission and transition probabilities
	needed for fitting the hidden Markov model.
}
\description{

	Computes emission probabilities and transition
	probabilities. See details.
}
\usage{
hmm.setup(object, states = paste("state", 1:length(copynumberStates),
	sep = ""), ICE = FALSE, copyNumber = TRUE, copynumberStates =
	0:4, EMIT.THR = -10, scaleSds = TRUE, verbose = TRUE,
	log.initial = log(rep(1/length(states), length(states))),
	normalIndex = 3, prGenotypeHomozygous = numeric(),
	prGenotypeMissing = rep(1/length(states), length(states)),
	pHetCalledHom = 0.001, pHetCalledHet = 0.995, pHomInNormal =
	0.8, pHomInRoh = 0.999, rohStates = logical(), trioHmm = FALSE, ...)
}

\arguments{
  \item{object}{
  The object must be one of the following classes inherited from
  \code{eSet} and defined in the R package oligoClasses:
  \code{SnpSet}, \code{oligoSnpSet}, \code{CopyNumberSet}, or
  \code{CNSet}.

}

  \item{states}{

  Vector of names for the hidden states.

}

  \item{ICE}{

  If the object is of class \code{oligoSnpSet} or \code{SnpSet} and
  the R package crlmm was used to call genotypes, the computed
  emission probabilities incorporate the confidence estimates of the
  genotype calls.

}

  \item{copyNumber}{
  
  Logical.  Whether to include information on copy number in the
  hidden Markov model.  If the object if of class \code{SnpSet}, this
  argument is set to FALSE.

}
  \item{copynumberStates}{
  
  Numerical vector with same length as the number of \code{states}.
  Each value corresponds to the latent copy number of the hidden
  state.  Note that copynumberStates must be specified on the
  appropriate scale.  If the copy number estimates have been
  log-transformed, the copynumberStates must be provided on the
  log-scale.  

}

  \item{EMIT.THR}{

  Single point outliers can cause the HMM to be jumpy. Emission
  probabilities below EMIT.THR are set to EMIT.THR.

}

  \item{scaleSds}{

  Logical.  For objects of class \code{CNSet}, sd estimates for total
  copy number are obtained by using \code{robustSds} function.

}
  \item{verbose}{

  Logical.  Verbose output during calculations.

}

  \item{log.initial}{

  Numeric vector of initial state probabilities (log-scale)
  corresponding to the hidden states.  Must be the same length as
  \code{states}.

}

  \item{normalIndex}{
  
  Integer.  

  \code{states[normalIndex]} should return the name of the 'normal'
  hidden state, which in general corresponds to copy number 2.  For
  instance, if the states were "hemizygousDeletion", "normal", and
  "amplification", normalIndex is 2.

}

 \item{prGenotypeHomozygous}{
  
  Numeric. The probability of a homozygous genotype call in each of
  the hidden states.  Ignored if ICE is TRUE.

}

  \item{prGenotypeMissing}{
  
  Numeric. The probability of a missing genotype for each hidden
  states.  
  
}

  \item{pHetCalledHom}{

  Numeric.  Probability of misclassifying a genotype call as
  heterozygous if the true genotype is homozygous.  Ignored unless ICE
  is TRUE.

}

  \item{pHetCalledHet}{

  Numeric.  Probability of correctly classifying a genotype call as
  heterozygous. Ignored unless ICE is TRUE.


}

  \item{pHomInNormal}{

  Numeric.  Probability of a homozygous genotype in a region without
  loss of heterozygosity.  Ignored unless ICE is TRUE.


}

  \item{pHomInRoh}{

  Numeric.  Probability of a homozygous genotype in a 'region of
  homozygosity'.  Ignored unless ICE is TRUE.


}

  \item{rohStates}{
  
  Logical vector.  TRUE corresponds to a hidden states in which
  regions of homozygosity are expected.  For instance, regions of
  homozygosity would be TRUE for hidden states corresonding to
  copy-neutral region of homozygosity (as my occur in a loss of
  heterozygosity region) and hemizygous deletions.

}

 \item{trioHmm}{

 Logical.  This option is experimental.  For Father-Mother-Offspring
 trios, we compute emission probabilities for biparental inheritance
 where the genotypes are informative. The hidden states correspond to
 biparental inheritance or non-biparental inheritance.  Regions of
 non-biparental inheritance can be used to quickly flag regions that
 are possibly de-novo deletions.

}
  \item{\dots}{
  Ignored.
}
}

\details{
	Details on the calculation of emission probabilities.
}

\value{

}
\author{
R. Scharpf
}

\seealso{
	\code{\link{robustSds}}
}
\examples{
}
\keyword{models}
\keyword{manip}

