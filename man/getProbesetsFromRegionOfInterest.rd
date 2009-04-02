
\name{getProbesetsFromRegionOfInterest}
\alias{getProbesetsFromRegionOfInterest}

\title{Get Probeset IDs From Region Of Interest}
\description{
Function that return the probesets located within a given region of the genome.
}
\usage{
    getProbesetsFromRegionOfInterest(annotation, chromosome, start, end, pythonPath=NULL, transcriptClustersFile=NULL, mpsToPsFile=NULL)
}
\arguments{
  \item{annotation}{A characther string giving the type of array for which probesets are needed. Will be used to load the .db file from bioconductor. Optional if mps to ps and transcript cluster file is given.}
  \item{chromosome}{The chromosome of interest. Should be given as a character string of the type "Chr1", "Chr2", "ChrY".}
  \item{start}{A character string with the start location of interest in the chromosome.}
  \item{end}{A character string with the end location of interest in the chromosome.}
  \item{pythonPath}{Optional character string with the path for Python software. This is only needed for exon arrays. If Python is in path this will be recognised automatically. Python can be downloaded from http://www.python.org.}
  \item{transcriptClustersFile}{The location of a transcript cluster file such as HuEx-1\_0-st-v2.na26.hg18.transcript.csv. These can be downloaded from http://www.affymetrix.com.}
  \item{mpsToPsFile}{The location of a transcript cluster file such as HuEx-1\_0-st-v2.r2.dt1.hg18.full.mps. These can be downloaded from http://www.affymetrix.com. Best to use the "full" type file.}
}
\value{A list of all probesets found in the given range on the given chromosome.}
\details{While working with 3'IVT type affymetrix arrays, for which .db files exist within the bioconductor environment, this 
function works in a pretty simple way. However if exon arrays are used, it will change to parsing the relevant library
files from the http://www.affymetrix.com website using a python script. Admittedly this is not optional from an R-only
point-of-view, but it works and its fast. The function will be updated when more R-centric ways of parsing exon array 
annotation data is available. 

Alternatively this data can just as well be retrieved from the web, but in some cases this function is faster an easier. 
}
\author{Lasse Folkersen}
\seealso{\code{\link{getLocalProbeIntensities}}, \code{\link{getMetaprobesetsFromRegionOfInterest}}}
\examples{

	\dontrun{
	#must supply mpsToPsFile and transcriptClustersFile for this to work
	probesets<-getProbesetsFromRegionOfInterest("notusedhere", chromosome=2, start="215889955", end="216106710", transcriptClustersFile=transcriptClustersFile, mpsToPsFile=mpsToPsFile)
	}
}
\keyword{documentation}
\keyword{utilities}

