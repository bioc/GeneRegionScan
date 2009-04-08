
\name{geneRegionScan}
\alias{geneRegionScan}
\alias{GeneRegionScan}
\alias{geneRegionScan,ExpressionSet-method}


\title{Gene Region Scan}
\description{
The top-level wrapper function that outputs as much data as possible, concerning one or a few genes. Refer to the functions
\link{plotCoexpression} and \link{plotOnGene} for further explanation of each part of the plot. 
}
\usage{
    geneRegionScan(object, gene, genomicData=NULL, probeData=NULL, label=NULL, genename=NULL, summaryType="median", yMax=NULL, testType=NULL,forcePValue=FALSE,verbose=TRUE,cutoff=0.2, directions="all",correlationCutoff=0.3, probeLevelInfo=c("probeid"))
    }
\arguments{
  \item{object}{A ProbeLevelSet object or a regular ExpressionSet object (in which case a probeData argument is required). See \link{getLocalProbeIntensities} and related functions on how to create a ProbeLevelSet.}
  \item{gene}{A number of gene sequences as DNAString, list of DNAString objects, character-vectors or \link[Biostrings]{readFASTA} outputs.}
  \item{genomicData}{Optional. If only one gene is specified this can be of the same form as given in \link{exonStructure}. If more than one gene is given, it must be a list, containing of one of these forms of the argument for each of the genes, in the same order.}
  \item{probeData}{Optional if a ProbeLevelSet is submitted as object argument. Otherwise, it must be a data frame with rownames corresponding to the featureNames of the ExpressionSet and a column named "sequence" with the probe sequences as character strings}
  \item{label}{An optional character string specifying a column name in the pData of the object. If this argument is given, the gene plot will be colour coded based on the different groups (factors) in the pData entry. If a summaryType other than 'dots' is selected the summarisation is done stratified by the different groups in the pData.}    
  \item{genename}{Optional character string specifying a gene name to include in the plot. If not included and a FASTA sequence is given, it will default to the name in the FASTA sequence. Otherwise it will default to 'Unknown genename'.}
  \item{summaryType}{Character string specifying one of the following summary methods: 'median', 'mean', 'quartiles' or 'dots' (i.e. no summary). Specifies how all the sample values or all the samples values in a group if 'label' is given, should be summarised. Defaults to 'median'.}
  \item{yMax}{Optional integer. If given, this value will be the maximal value on the y-axis. This is useful if a few outlier probes have very high intensity values, as the default is to set the yMax to the maximal intensity value.}
  \item{testType}{Optional character string, defining a statistic procedure to identify especially interesting probes. Can be either 'linear model', 'students' or 'wilcoxons'. If given, a label must also be specified. In this case the \link{plotStatistics} function will be called and probes that are significantly changed between the groups in label at the P-value set in cutoff (see cutoff argument) will be circled.}
  \item{forcePValue}{Logical. Is used if the testType argument is used. If TRUE all significantly changed probes have P-value given on the plot. If FALSE, only plots with less than 10 significant probes write P-values. Plots can become very cluttered with data if set to TRUE}
  \item{verbose}{TRUE or FALSE}
  \item{cutoff}{Integer specifying at what p-value probes should be circled when using the 'testType' variable. Defaults to 0.2. For cutoffs higher than 0.05, all probes with P >0.05 will be circled in grey instead of black.}
  \item{directions}{A character vector of the matching-directions that should be scanned (which combinations of complementary and reverse). Defaults to "all" which is shorthand for all possible directions, but can take anything from: c("matchForwardSense", "matchForwardAntisense", "matchReverseSense", "matchReverseAntisense")}
  \item{correlationCutoff}{A number between 0 and 1. The limit at which Pearson correlation (in absolute values) should not be plotted below. Defaults to 0.5}
  \item{probeLevelInfo}{The information about each probe to include in the plot. Should be a vector of one or more of the following elements: probeid, probesetid, sequence. Default is only probeid.}
  
    
}
\value{No value, but plots the local expression levels of each probe found in the submitted sequence of the gene or 
genes as a function of its location along this sequence on the top half of a pdf-file. On the bottom half of the pdf-file 
it will plot the pairwise correlations between all probes found in the sequences. The pdf file will be named "report\_" + title
of ExpressionSet or ProbeLevelSet}
\details{This function is a wrapper around \link{plotOnGene} and \link{plotCoexpression}. The output of plotOnGene is 
included in the top half of the plot and the output of plotCoexpression will be included in the bottom half of the plot.
Refer to each of these functions for more detailed help. 

In general, this function gives an overview of how intensity values of individual probes on a microarray are in relation
to an actual gene or set of genes in a region of the genome. 

The function has only been tested with up to four genes at the same time. A plot with more genes would probably also be 
too complicated to interpret with this method. In addition, the alignment of the top and bottom plots also becomes somewhat
difficult with more genes. This alignment is also the reason why the function can not export to any device.

See \link{getLocalProbeIntensities} for more info on how to obtain ProbeLevelSets.
}
\author{Lasse Folkersen}
\seealso{\code{\link{plotCoexpression}}, \code{\link{plotOnGene}}, \code{\link{getLocalProbeIntensities}}}
\examples{

	data(exampleProbeLevelSet)
	
	#simple:
	geneRegionScan(exampleProbeLevelSet,mrna)
	
	#more complicated - note that we slice the mrna to simulate comparing two different isoforms
	gene1<-DNAString(mrna[[1]]$seq)[1:1000]
	gene2<-DNAString(mrna[[1]]$seq)[1500:3000]
	
	geneRegionScan(exampleProbeLevelSet, list(gene1,gene2), genomicData=list(genomic,genomic), label="genotype3", summaryType="mean",
    testType="linear model", forcePValue=TRUE, cutoff=0.1, directions="all", correlationCutoff=0.6,
    probeLevelInfo=c("probeid","sequence"))
}
\keyword{documentation}
\keyword{utilities}

