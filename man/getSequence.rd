\name{getSequence}
\alias{getSequence}
\alias{getSequence,ProbeLevelSet-method}

\title{Get Sequence from a ProbeLevelSet}
\description{
Function to retrieve the sequences of feature in a ProbeLevelSet.
}
\usage{
    getSequence(object, id)
}
\arguments{
  \item{object}{A ProbeLevelSet object}
  \item{id}{Optional character vector with the featureNames of the sequences of interest}
}
\value{A character vector with the sequences of all probes in the ProbeLevelSet. The names of the vector are set to match the sequences.}

\author{Diego Diez}
\seealso{\code{\link{ProbeLevelSet}}}
\examples{

	data(exampleProbeLevelSet)
	getSequence(exampleProbeLevelSet)
	
}
\keyword{documentation}
\keyword{utilities}

