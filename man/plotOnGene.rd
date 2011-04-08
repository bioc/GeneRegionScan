
\name{plotOnGene}
\alias{plotOnGene}
\alias{plotOnGene,ExpressionSet-method}

\title{Plot probe level data on a gene}
\description{
Function that will investigate the probe level intensity of probes as a function of their location in a gene.
}
\usage{
    plotOnGene(object, gene, probeData=NULL, label=NULL, genename=NULL, summaryType="median",
    interval=NULL, ylim=NULL, testType=NULL, forcePValue=FALSE, verbose=TRUE, cutoff=0.2, directions="all", ylab="expression")
}
\arguments{
  \item{object}{A ProbeLevelSet object or a regular ExpressionSet object (in which case a probeData argument is required). See \link{getLocalProbeIntensities} and related functions on how to create a ProbeLevelSet.}
  \item{gene}{A number of gene sequences as DNAstring, vectors of DNAStrings, character-vectors or \link[Biostrings]{readFASTA} outputs.}
  \item{probeData}{Optional if a ProbeLevelSet is submitted as object argument. Otherwise, it must be a data frame with rownames corresponding to the featureNames of the ExpressionSet and a column named "sequence" with the probe sequences as character strings}
  \item{label}{An optional character string specifying a column name in the pData of the object. If this argument is given, the gene plot will be colour coded based on the different groups (factors) in the pData entry. If a summaryType other than 'dots' is selected the summarisation is done stratified by the different groups in the pData. It can be a numeric or integer entry, but it will be coerced to factors.}    
  \item{genename}{Optional character string specifying a gene name to include in the plot. If not included and a FASTA sequence is given, it will default to the name in the FASTA sequence. Otherwise, it will default to 'Unknown genename'.}
  \item{summaryType}{Character string specifying one of the following summary methods: 'median', 'mean', 'quartiles' or 'dots' (i.e. no summary). Specifies how all the sample values or all the samples values in a group if 'label' is given, should be summarised. Defaults to 'median'.}
  \item{interval}{Optional vector of two integers of bp positions. If given, the plot will only include the sequence from gene in the given interval. The x-axis annotation is preserved from original, so this is useful for zooming on specific regions.}
  \item{ylim}{Optional two integers. If given, this value will be the minimal and maximal value on the y-axis. This is useful if a few outlier probes have very high intensity values, as the default is to set the ylim from the maximal intensity value.}
  \item{testType}{Optional character string, defining a statistic procedure to identify especially interesting probes. Can be either 'linear model', 'students' or 'wilcoxons'. If given, a label must also be specified. In this case the \link{plotStatistics} function will be called and probes that are significantly changed between the groups in label at the P-value set in cutoff (see cutoff argument) will be circled.}
  \item{forcePValue}{Logical. Is used if the testType argument is used. If TRUE all significantly changed probes have P-value given on the plot. If FALSE, only plots with less than 10 significant probes write P-values. Plots can become very cluttered with data if set to TRUE}
  \item{verbose}{TRUE or FALSE}
  \item{cutoff}{Integer specifying at what p-value probes should be circled when using the 'testType' variable. Defaults to 0.2. For cutoffs higher than 0.05, all probes with P >0.05 will be circled in grey instead of black.}
  \item{directions}{A character vector of the matching-directions that should be scanned (which combinations of complementary and reverse). Defaults to "all" which is shorthand for all possible directions, but can take anything from: c("matchForwardSense", "matchForwardAntisense", "matchReverseSense", "matchReverseAntisense")}
  \item{ylab}{Label of Y-axis as in default plots}
  
    
}
\value{No value, but plots the local expression levels relations of each probe found in the submitted gene sequence as a function of its location along this sequence. Various statistics and summarizations on pdata can be employed, as specified in details.}
\details{At the very least, this function takes a ProbeLevelSet or an ExpressionSet + probeData and the sequence of a gene. It then
compares the probe sequences given in the ProbeLevelSet or the probeData variable with the sequence of the gene given. Any probes with
sequences found in the gene will be plotted, with their intensity level on the y-axis and their location in the gene on the x-axis.
 If no further arguments are given this gives a view of relative expression levels along the length of the gene, and can be used
 to investigate which exons are actively transcribed in the sample and which are not. An important argument that can be used for further
 investigation is the 'label' argument which specifies a column in the pData of the ExpressionSet. In this case the plots will be
 stratified by the factors specified in this column (so giving labels with numerical or Date class data will not work). This can be
 very useful when investigating how different sample conditions affect various regions of a gene. A transcript isoform that is relatively
 upregulated in a diseased state will for example not be discovered if a probeset or metaprobeset covering the entire gene is used
 to summarize the data, since the average expression intensity for the gene will remain constant. Using the plotOnGene function, however,
 and specifying a case / control label will reveal tendencies for probes at certain exon locations to have relation to this label.
 The testType argument further supports this functionality by providing statistical testing and highlighting of probes that correlate
 significantly to the given label. In the case / control example, a student's t-test would highlight all probes that matched with the
 exons of the gene that was only found in the disease-specific transcript isoform. When interpreting the data it is suggested that
 specific attention is paid to the pattern of probes in the same exon. A single probe with a P-value < 0.05 might be a false positive
 caused by chance or by cross hybridization of the probe sequence to something else. A range of probes in the same exon that all show
 P-values below or close to 0.05, however, is much more likely to be an actual case of a transcript isoform having this particular 
 exon or exons being regulated between the groups in the label. Exon structure can be easily plotted on the graph using the
 \link{exonStructure} function.
 
 A special case is the search for SNPs which have effect on expression levels or variable splicing. The testType argument
 'linear model' is designed for this. The linear model calls the internal function \link{doProbeLinear} which assign each of the levels
 in the 'label' column of the pdata a value between 1 and the number of levels, in the order in which they are sorted. For genotypes
 given as "AA", "AB", "BB" character strings this will give "AA" = 1, "AB" = 2 and "BB" = 3. The doProbeLinear then calculates a
 linear model between the intensity values and these numbers, and returns the P-value. In this case, low P-values can be interpreted
 as a case where the heterozygote samples have intermediary expression levels between the two homozygotes. This is the case that
 can be expected to be seen if the nucleotide type of SNP does in fact have any effect on the  mRNA concentration levels in the 
 sample.
}
\author{Lasse Folkersen}
\seealso{\code{\link{geneRegionScan}}, \code{\link{plotCoexpression}}}
\examples{

	data(exampleProbeLevelSet)
	
	plotOnGene(exampleProbeLevelSet, mrna, summaryType="dots", interval=c(500,1000))
	
	plotOnGene(exampleProbeLevelSet, mrna, label="genotype3", testType="linear model")
	exonStructure(mrna, genomic)
}
\keyword{documentation}
\keyword{utilities}

