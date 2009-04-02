
\name{checkForFileInPath}
\alias{checkForFileInPath}

\title{Check For File In Path}
\description{
Internal function used to check for files in path.
}
\usage{
    checkForFileInPath(filenames)
}
\arguments{
  \item{filenames}{A character vector of filenames to check for - c("xxxx.exe","xxxx") for example}
}
\value{A character vector with the expanded paths and filenames to the found files.}
\details{
This function is not meant to be called by the user so it is not in the namespace.
} 
 
\author{Lasse Folkersen}
\seealso{ \link{plotOnGene} }
\examples{

	
	GeneRegionScan:::checkForFileInPath(c("python","python.exe"))
	
}
\keyword{documentation}
\keyword{utilities}

