\name{locusLevelData}
\alias{locusLevelData}
\docType{data}
\title{Basic data elements required for the HMM}
\description{
  This object is a list containing the basic data elements required for the HMM
}
\usage{data(locusLevelData)}
\format{
  A list
}
\details{
  
  The basic assay data elements that can be used for fitting the HMM
  are:

  1.  a mapping of platform identifiers to chromosome and physical
  position

  2. (optional) a matrix of copy number estimates

  3. (optional) a matrix of confidence scores for the copy number
  estimates (e.g., inverse standard deviations)
  
  4. (optional) a matrix of genotype calls

  5. (optional) CRLMM confidence scores for the genotype calls

  At least (2) or (4) is required.  The locusLevelData is a list that
  contains (1), (2), (4), and (5).
  
}
\source{
  
  A HapMap sample on the Affymetrix 50k platform.  Chromosomal
  alterations were simulated.  The last 100 SNPs on chromosome 2 are, in
  fact, a repeat of the first 100 SNPs on chromosome 1 -- this was added
  for internal use.
  
}

\examples{
data(locusLevelData)
str(locusLevelData)
}
\keyword{datasets}
