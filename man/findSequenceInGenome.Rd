
\name{findSequenceInGenome}
\alias{findSequenceInGenome}

\title{Find a sequence in genome}
\description{
Wrapper around \link[Biostrings:matchPDict-exact]{matchPDict} that will
accept a list of sequences and check if they are present in a given
genome. Takes a long time to run. 
}
\usage{
    findSequenceInGenome(sequences, genome="BSgenome.Hsapiens.UCSC.hg18", verbose=TRUE,
    directions=c("matchForwardSense", "matchForwardAntisense", "matchReverseSense", "matchReverseAntisense"))
}
\arguments{
  \item{sequences}{vector of character strings to scan. Should only contain A, C, G and T. Will be converted to DNAString.}
  \item{genome}{character string with the name of the BSGenome in which sequences should be found. Defaults to the human genome.}
  \item{verbose}{TRUE or FALSE.}
  \item{directions}{character string with elements from c("matchForwardSense", "matchForwardAntisense", "matchReverseSense", "matchReverseAntisense"). Defines which directions (complementary and reverse mirrorings) that should be scanned. Defaults to all directions.}
}
\value{A data frame with a row for each identified match. Columns are "entrynumber","hitposInChr","chr", and "sequence" describing, respectively: the index of the sequence match, the position in the chromosome at which it was found, which chromosome it was found on, the sequence itself}
\details{This function will take quite a while to run, so if you have a many sequences, overnight runs are recommended. BSgenome contains some alternative versions of chromosomes. They are marked with an underscore. This function automatically disregards chromosome names with an underscore, and this is known to work for the human genome. Nevertheless, check the output printed to terminal if all chromosomes are included.}
\author{Lasse Folkersen}
\seealso{\code{\link[BSgenome:BSgenome-class]{BSgenome}},
  \code{\link[Biostrings]{matchPDict}},
  \code{\link{excludeDoubleMatchingProbes}}} 
\examples{



	\dontrun{
	#you can run this, but it takes quite a lot of time
	example<-findSequenceInGenome("CTGGCGAGCAGCGAATAATGGTTT")
	}
}
\keyword{documentation}
\keyword{utilities}

