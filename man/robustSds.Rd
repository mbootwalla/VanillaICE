\name{robustSds}
\alias{robustSds}
\title{
	Calculate robust estimates of the standard deviation
}
\description{
	Uses the median absolute deviation (MAD) to calculate robust estimates of the standard deviation 
}
\usage{
robustSds(x, takeLog = FALSE, ...)
}
\arguments{
  \item{x}{
	A matrix of copy number estimates. Rows are features, columns are samples.
}
  \item{takeLog}{
  Whether to log-transform the copy number estimates before computing
  robust sds
}
  \item{\dots}{
  additional arguments to \code{rowMedians}
}
}
\details{
	
	For matrices \code{x} with 4 or more samples, the row-wise MAD
	(SNP-specific sds) are scaled by sample MAD / median(sample
	MAD).

	If the matrix has 3 or fewer samples, the MAD of the sample(s)
	is returned.

}
\value{
	Matrix of standard deviations.
}
\examples{
	data(locusLevelData)
	sds <- robustSds(locusLevelData[["copynumber"]]/100,
			 takeLog=TRUE)


}
\keyword{ manip }

