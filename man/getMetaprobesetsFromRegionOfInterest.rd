
\name{getMetaprobesetsFromRegionOfInterest}
\alias{getMetaprobesetsFromRegionOfInterest}

\title{Get Meta Probe Set IDs From Region Of Interest}
\description{
Function that return the meta probe set ids located within a given region of the genome.
}
\usage{
	getMetaprobesetsFromRegionOfInterest(annotation, chromosome, start, end, pythonPath=NULL, transcriptClustersFile=NULL)
}
\arguments{
  \item{annotation}{A characther string giving the type of array for which probesets are needed. This is used to save the location of the transcriptClustersFile.}
  \item{chromosome}{The chromosome of interest. Should be given as a character string of the type "Chr1", "Chr2", "ChrY".}
  \item{start}{A character string with the start location of interest in the chromosome.}
  \item{end}{A character string with the end location of interest in the chromosome.}
  \item{pythonPath}{Optional character string with the path for Python software. This is only needed for exon arrays. If Python is in path this will be recognised automatically. Python can be downloaded from http://www.python.org.}
  \item{transcriptClustersFile}{The location of a transcript cluster file such as HuEx-1\_0-st-v2.na26.hg18.transcript.csv. These can be downloaded from http://www.affymetrix.com.}
}  
\value{A list of all probesets found in the given range on the given chromosome.}
\details{Since meta probe sets are only found in exon arrays, this function does not work with 3' IVT array files.
It parses the location files from the http://www.affymetrix.com website using a python script. Admittedly this is not optional from an R-only
point-of-view, but it works and its fast. The function will be updated when more R-centric ways of parsing exon array 
annotation data is available. 

Alternatively this data can just as well be retrieved from the web, but in some cases this function is faster an easier. 
}
\author{Lasse Folkersen}
\seealso{\code{\link{getLocalProbeIntensities}}, \code{\link{getProbesetsFromRegionOfInterest}}}
\examples{

	\dontrun{
	#must supply transcriptClustersFile for this to work
	metaprobesets<-getMetaprobesetsFromRegionOfInterest("notusedhere", chromosome=2, start="215889955", end="216106710", transcriptClustersFile=transcriptClustersFile)
	}
}
\keyword{documentation}
\keyword{utilities}

