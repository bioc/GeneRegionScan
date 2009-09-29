
\name{excludeDoubleMatchingProbes}
\alias{excludeDoubleMatchingProbes}
\alias{excludeDoubleMatchingProbes,ProbeLevelSet-method}


\title{Exclude ProbeLevelSet probes that match more than once in genome}
\description{
Function that will remove probes from ProbeLevelSet if they have more than one match in a given genome. 
}
\usage{
    excludeDoubleMatchingProbes(object, genome="BSgenome.Hsapiens.UCSC.hg18", verbose=TRUE,
    directions=c("matchForwardSense", "matchForwardAntisense", "matchReverseSense", "matchReverseAntisense"),
    previousData = NULL)
}
\arguments{
  \item{object}{A ProbeLevelSet class object.}
  \item{genome}{character string with the name of the BSGenome in which sequences should be found. Defaults to the human genome.}
  \item{verbose}{TRUE or FALSE.}
  \item{directions}{character string with elements from c("matchForwardSense", "matchForwardAntisense", "matchReverseSense", "matchReverseAntisense"). Defines which directions (complementary and reverse mirrorings) that should be scanned. Defaults to all directions.}
  \item{previousData}{Optional: The output from a call to findSequenceInGenome. If given the scanning will be skipped, and the probes will be omitted directly. Useful in cases were datasets from the same region needs to be processed.}
}
\value{The ProbeLevelSet class object provided as argument, with all double matching probes removed. Double matching probes are probes whose sequence are found twice or more in the genome. In addition, the output of the matching investigation is saved in the notes of the ProbeLevelSet and can be further examined for information on the locations of the probe sequences in relation to the BSgenome sequences.}
\details{This function will take quite a while to run, so if you have many sequences, overnight runs are recommended. 
BSgenome contains some alternative versions of chromosomes. They are marked with an underscore. This function 
automatically disregards chromosome names with an underscore, and this is known to work for the human genome. Nevertheless, 
check the output printed to terminal if all chromosomes are included. The function is a wrapper around findSequenceInGenome, 
which can be used for purposes that are more flexible (although that function really is just following the example in the BSgenome package)

At present, there is no functionality to check matches with known SNP or known splice forms taken into account.
}
\author{Lasse Folkersen}
\seealso{\code{\link{findSequenceInGenome}},
  \code{\link[BSgenome:BSgenome-class]{BSgenome}},
  \code{\link{excludeDoubleMatchingProbes}}} 
\examples{

	\dontrun{
	#you can run this, but it takes a lot of time
	probelevelsetwithnodoubles<-excludeDoubleMatchingProbes(probelevelset)
	}
}
\keyword{documentation}
\keyword{utilities}

