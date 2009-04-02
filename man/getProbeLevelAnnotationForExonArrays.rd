
\name{getProbeLevelAnnotationForExonArrays}
\alias{getProbeLevelAnnotationForExonArrays}

\title{Get ProbeLevel Annotation for Exon Arrays}
\description{
Internal function that will return exon probe sequence and probe id for all probes in a given list of probe sets.
}
\usage{
    getProbeLevelAnnotationForExonArrays(vectorOfProbesets, pgfPath)
}
\arguments{
  \item{vectorOfProbesets}{A character string with probeset IDs.}
  \item{pgfPath}{The path of a pgf file for the exon array type of interest. Can be downloaded from www.affymetrix.com.}
}
\value{A dataframe with a row for each probe in the submitted probeset. The columns with the names "probeset\_name" and "sequence" contains this.}
\details{This function makes a call to readPgf in the affxparser package to get the sequence information. The readPgf will then load the entire pgf file
and that might be quite memory intensive. However it is not possible to extract the indices of the probeset names without doing this. The function is primarily
intended to be called by \link{getLocalProbeIntensities}.

The call to \link[affxparser]{readPgf} sometimes gives a is.na() warning. The reason for this is not known, but it does not seem
to affect performance.
}
\author{Lasse Folkersen}
\seealso{\code{\link[affxparser]{readPgf}}, \code{\link{getLocalProbeIntensities}}, \code{\link{plotOnGene}}}
\examples{

	\dontrun{
	#must supply pgf file for this to work
	getProbeLevelAnnotationForExonArrays("43254543", pgfPath="~/somewhere/some.pgf")
	}
}
\keyword{documentation}
\keyword{utilities}

