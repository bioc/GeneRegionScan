\name{genomic}
\docType{data}
\alias{genomic}
\title{genomic sequence of gene in ProbeLevelSet}
\description{
This is the DNA sequence for fibronectin 1 isoform 2 preproprotein, which is used in the exampleProbeLevelSet.
It has been split up so each exon in the gene is described in one FASTA entry.


  }
\usage{data(exampleProbeLevelSet)}
\format{
	It has been read from FASTA format using the \link[Biostrings:readDNAStringSet]{readDNAStringSet} function.
	The FASTA was downloaded from http://genome.ucsc.edu, specifying "One FASTA record per region"
	in the "Sequence Retrieval Region Options".
}
\author{Lasse Folkersen}
\keyword{datasets}