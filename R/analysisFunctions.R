
setGeneric("doProbeLinear", function(object,label,testType="linear model") standardGeneric("doProbeLinear"))
setMethod("doProbeLinear", "ExpressionSet",
		function(object,label,testType="linear model"){
			# Function to calculate if genes are interesting using linear models 
			# takes an expressionset and a label
			# outputs a matrix with the following data for each probe set "p-value","slope","intercept","slope/mean", and "anova"
			# This function is not meant to be run by the user. Only to be called from elsewhere.
			
			oldc <- Sys.getlocale("LC_COLLATE")
			on.exit(Sys.setlocale("LC_COLLATE", oldc))
			Sys.setlocale("LC_COLLATE", "C")
			
			expressionset <- object
			
			if(!testType == "linear model"){stop(paste("The testType",testType,"is not recognised"))}
			
			factors <- sort(levels(pData(expressionset)[,label]))
			
			genotype <- vector(mode="integer",length=ncol(expressionset))
			message_string <- ""
			for(i in 1:length(factors)){
				genotype[pData(expressionset)[,label] == factors[i]] <- i
				message_string <- paste(message_string,factors[i],"=",i) 
			}
			
			
			
			
			
			expressionset <- expressionset[,genotype != 0] # Remove arrays where the label does not exists. Next four lines re-calculate the genotype vector.
			genotype <- vector(mode="integer",length=ncol(expressionset))
			for(i in 1:length(factors)){
				genotype[pData(expressionset)[,label] == factors[i]] <- i	
			}
			
			correlation_list <- matrix(ncol=4,nrow=nrow(expressionset))
			colnames(correlation_list) <- c("p-value","slope","intercept","slope/mean")
			rownames(correlation_list) <- featureNames(expressionset)
			for(i in 1:nrow(expressionset)){
				gene_expression <- exprs(expressionset)[i,]
				model <- lm(gene_expression~genotype)
				correlation_list[i,] <- c(
						summary(model)[["coefficients"]]["genotype","Pr(>|t|)"],
						summary(model)[["coefficients"]]["genotype","Estimate"],
						summary(model)[["coefficients"]]["(Intercept)","Estimate"],
						summary(model)[["coefficients"]]["genotype","Estimate"]/mean(exprs(expressionset)[i,])
				)
			}
			return(correlation_list)
		})


setGeneric("doProbeTTest", function(object,label,testType="students") standardGeneric("doProbeTTest"))
setMethod("doProbeTTest", "ExpressionSet",
		function(object,label,testType="students"){
			#Function to calculate if genes are interesting using students and wilcoxons
			#takes an expressionset and a label
			#testType
			#	options are:
			#		"students" - parametric students t.test
			#		"wilcoxon" - non-parametrix wilcoxon test
			#outputs a matrix with the following data for each probe set "p-value","fold-change"
			expressionset <- object
			
			factors <- levels(pData(expressionset)[,label])
			
			
			
			if(!(testType == "students"|testType == "wilcoxon")){stop(paste("The testType",testType,"is not recognised"))}
			if(length(factors) != 2){stop(paste(length(factors),"different factors found in",label,"- the algorithm has not been designed for anything else than 2 factors"))}
			
			
			group_one <- which(pData(expressionset)[,label] == factors[1])
			group_two <- which(pData(expressionset)[,label] == factors[2])
			
			if(length(group_one)<2 | length(group_two)<2)stop("One of the groups compared had less than two members. A T-test can't be calculated under these conditions")
			
			correlation_list <- matrix(ncol=2,nrow=nrow(expressionset))
			colnames(correlation_list) <- c("p-value","fold-change")
			rownames(correlation_list) <- featureNames(expressionset)
			for(i in 1:nrow(expressionset)){
				
				if(testType == "students")p_value <- t.test(exprs(expressionset)[i,group_one],exprs(expressionset)[i,group_two],na.rm=TRUE)[["p.value"]]
				if(testType == "wilcoxon")p_value <- wilcox.test(exprs(expressionset)[i,group_one],exprs(expressionset)[i,group_two],na.rm=TRUE)[["p.value"]]
				fold_change <- mean(exprs(expressionset)[i,group_one],na.rm=TRUE)/mean(exprs(expressionset)[i,group_two],,na.rm=TRUE)
				if(fold_change<1){
					fold_change <- -(1/fold_change)
				}
				correlation_list[i,] <- c(p_value,fold_change)
			}
			return(correlation_list)
		})

setGeneric("plotStatistics", function(object,probeData,label,summaryType,testType,interval=NULL,forcePValue=FALSE,verbose=TRUE,positionVector,ylim,xlim,cutoff=0.2) standardGeneric("plotStatistics"))
setMethod("plotStatistics", "ExpressionSet",
		function(object,probeData,label,summaryType,testType,interval=NULL,forcePValue=FALSE,verbose=TRUE,positionVector,ylim,xlim,cutoff=0.2){
			# Function to plot statistics
			# Not meant to be called by users. Only other functions
			# will take decorate a plot of probe level expression data with markings of particularly significant relations.
			# The cutoff will determine at which P-values probes are highlighted. P-values higher than 0.05 are always highlighted grey instead of black.
			# Other arguments are explained in plotOnGene
			expressionset <- object
			if(missing(label))stop("Statistical testing can not be done if no label is given") 
			if(!label %in% colnames(pData(expressionset)))stop(paste("The label",label,"is not found in the expressionset"))
			
			factors <- levels(pData(expressionset)[,label])
			if(length(factors) == 1)stop(paste("The label",label,"only contains one type of variable, so no statistics can be done"))
			if(length(factors) > 3)stop(paste("The label",label,"contains more than 3 types of variable, no statistic tests have been implemented to analyse that"))
			
			if(missing(summaryType))stop("summaryType is missing from the plotStatistics function so it is not possible to calculate where to highlight probes")
			
			if(missing(testType))stop("A testType argument is needed to know which statistical test to calculate")
			if(testType == "student")testType <- "students"
			if(testType %in% c("linear","linear_model","linear models","linear models"))testType <- "linear model"
			if(testType == "wilcoxons")testType <- "wilcoxon"
			
			interval_length <- interval[2]-interval[1]
			
			if(length(factors) == 2|length(factors) == 3){
				if(length(factors) == 2){
					if(!((testType == "students")|(testType == "wilcoxon"))){
						warning(paste("The label given codes two variables:",paste(factors,collapse=", ")," - only the students and wilcoxon test types can be used for that. Test type:",testType,"was not recognised as one of them. Student's t-test has been used per default"))
						testType <- "students"
					}
					
					correlation_list <- doProbeTTest(expressionset[names(positionVector),],label=label,testType=testType)
				}
				if(length(factors) == 3){
					if(testType != "linear model"){
						warning(paste("The label given codes three variables:",paste(factors,collapse=", ")," - only the linear model test type can be used for that. Test type:",testType,"was not recognised as this. The linear model has been used per default"))
						testType <- "linear model"
					}
					correlation_list <- doProbeLinear(expressionset[names(positionVector),],label=label,testType=testType)
				}
				
				mtext(paste("Probes that contain significant correlations to variable according to a",testType,"test are circled"),padj=3.4,cex=0.7)
				
				significant_probes <- rownames(subset(correlation_list,correlation_list[,"p-value"]<cutoff))
				
				
				if(length(significant_probes) > 0){
					for(significant_probe in significant_probes){
						y <- vector()
						
						x <- positionVector[significant_probe]+interval[1]	
						
						for(j in 1:length(factors)){
							expression_vector <- exprs(expressionset[significant_probe,])[pData(expressionset)[,label] == factors[j]]
							if(summaryType == "median"){y[j] <- median(expression_vector,na.rm=TRUE)}
							if(summaryType == "mean"){y[j] <- mean(expression_vector,na.rm=TRUE)}
							if(summaryType == "quartiles"){
								y[j] <- quantile(expression_vector,na.rm=TRUE)[2]
								y[j+length(factors)] <- quantile(expression_vector,na.rm=TRUE)[3]
							}
							if(summaryType == "dots"){
								y[j] <- max(expression_vector,na.rm=TRUE)
								y[j+length(factors)] <- min(expression_vector,na.rm=TRUE)
							}
							p_value <- correlation_list[significant_probe,"p-value"]
						}
						
						ellipse_center <- c(x,min(y)+(max(y)-min(y))/2)
						ellipse_width <- (interval_length/40)
						ellipse_height <- (max(y)-min(y))*1.5
						
						yheight <- ylim[2]-ylim[1]
						xwidth <- xlim[2]-xlim[1]
						
						if(ellipse_height*xwidth<ellipse_width*yheight) ellipse_height <- ellipse_width*(yheight/xwidth)
						
						ellipse <- matrix(ncol=2,nrow=100)
						for(l in 1:100){
							radians <- (l/100)*2*pi
							ellipse[l,1] <- (cos(radians)*ellipse_width*0.5)+ellipse_center[1]
							ellipse[l,2] <- (sin(radians)*ellipse_height*0.5)+ellipse_center[2]
							
						}
						
						if(p_value<0.05)lines(ellipse,col="black")
						if(p_value >= 0.05)lines(ellipse,col="grey")
						
						#points(x=x,y=c,cex=size)
						
						if(verbose)print(paste("Significant correlation found at",x,"having a p-value of:",signif(p_value,digits=3)," - this is probe",significant_probe))
						#if(verbose)writeLines(paste("Significant correlation found at",x,"having a p-value of:",signif(p_value,digits=3)," - this is probe",significant_probe),con=stats_file)
						
						if(length(significant_probes)<10|forcePValue){
							if(p_value<0.05)text(x=ellipse_center[1]+2*ellipse_width,y=ellipse_center[2]+0.5*ellipse_height,label=paste("P =",signif(p_value,digits=2)),cex=0.6,col="black")
							if(p_value >= 0.05)text(x=ellipse_center[1]+2*ellipse_width,y=ellipse_center[2]+0.5*ellipse_height,label=paste("P =",signif(p_value,digits=2)),cex=0.6,col="grey")
						}
					}
				}
			}
		})




setGeneric("findProbePositions", function(object, gene, probeData=NULL,interval=NULL,directions="all",verbose=TRUE) standardGeneric("findProbePositions"))
setMethod("findProbePositions", "ExpressionSet",
		function(object, gene, probeData=NULL,interval=NULL,directions="all",verbose=TRUE){
			# This function figures out the position of probes based on their sequence and an input gene sequence
			#  * "expressionset"/"object" the expression of all probes and their associated pData. 
			#  * "gene" in the form of a loaded FASTA sequence, character-sequence or DNAString instance
			#  * "probeData" (Optional if ProbeLevelSet is given) a description of all features. It requires a data frame with the featureNames as rownames and a column "sequence" containing their sequence.
			#
			# Optional switches are:
			# * interval: Optional interval of the gene to zoom on. If numbers outside of the gene length is given, they're truncated to end with the gene
			# * verbose: Print a lot of outputtext?
			# * directions: A character vector of the matching-directions that should be scanned. Defaults to "all", which is shorthand for all the remaining directions: c("matchForwardSense","matchForwardAntisense","matchReverseSense","matchReverseAntisense")
			
			
			
			expressionset<-object
			if(!all(directions %in% c("matchForwardSense","matchForwardAntisense","matchReverseSense","matchReverseAntisense","all")))stop(paste("The directions variable:",paste(directions,collapse=", "),"was not recognised. The directions must be in \"matchForwardSense\",\"matchForwardAntisense\",\"matchReverseSense\",\"matchReverseAntisense\",\"all\""))
			
			if("all" %in% directions){
				directions <- c("matchForwardSense","matchForwardAntisense","matchReverseSense","matchReverseAntisense")
			}
			
			
			
			#This part choses if probe data is from the pData(featureData(x)) or from user-defined probe data (ie. a data frame with sequences)	
			if(is.null(probeData)){
				if(ncol(pData(featureData(expressionset))) == 0){
					stop(paste("No featuredata found in",expressionset@experimentData@title,"and no probe data explicitly given"))
				}else{
					probeData <- pData(featureData(expressionset))
					if(verbose)print(paste("Probe data for",nrow(probeData),"probes found as featureData in the expressionset",expressionset@experimentData@title))
				}
			}
			if(class(probeData[,"sequence"]) != "character")stop(paste("The 'sequence' column of the probeData had data of class",class(probeData[,"sequence"]),"and not character as required"))
			
			if(length(rownames(probeData)) == length(featureNames(expressionset))){
				if(!all(rownames(probeData) %in% featureNames(expressionset)))warning(paste("Not all probe ids given in the probe_level_annotation data in the notes section() of",expressionset@experimentData@title,"were found in the featureNames (probe ids) in the expressionset"))
				if(!all(featureNames(expressionset) %in% rownames(probeData)))warning(paste("Not all featureNames (probe ids) in the expressionset,",expressionset@experimentData@title,"were found in the probe ids given in the probe_level_annotation data in the notes section()"))
			}else{
				warning(paste("The expressionset,",expressionset@experimentData@title,"had",length(featureNames(expressionset)),"features in the exprs-part, but",length(rownames(probeData)),"entries in the probe_level_annotation of the notes() section. This might indicate corrupted data"))
			}
			
			
			gene <- readGeneInput(gene)
			if(length(gene) > 1){
				warning(paste("The input specified",length(gene),"genes, but the function only accepts one. The first one,",gene[[1]][["desc"]],"was used"))
				gene <- list(gene[[1]])
			}
			genename <- gene[[1]][["desc"]]
			sequence <- DNAString(gene[[1]][["seq"]])
			
			if(verbose)print(paste("Investigating",length(sequence),"bp sequence"))
			
			probe <- DNAStringSet(probeData[,"sequence"])
			
			
			positionVector <- vector()
			
			if("matchForwardSense" %in% directions){
				probe_pdict <- PDict(probe)
				matchForwardSense <- matchPDict(probe_pdict,sequence)
			}
			if("matchForwardAntisense" %in% directions){
				probe_pdict <- PDict(complement(probe))
				matchForwardAntisense <- matchPDict(probe_pdict,sequence)
			}
			if("matchReverseSense" %in% directions){ 
				probe_pdict <- PDict(reverse(probe))
				matchReverseSense <- matchPDict(probe_pdict,sequence)
			}
			if("matchReverseAntisense" %in% directions){
				probe_pdict <- PDict(reverseComplement(probe))
				matchReverseAntisense <- matchPDict(probe_pdict,sequence)
			}
			
			for(match_direction_name in directions){
				match_direction <- get(match_direction_name)
				if(sum(elementLengths(match_direction) > 1) > 0){
					stop(paste(sum(elementLengths(match_direction) > 1),"of the probes matched more than one time in the sequence"))
				}
				if(sum(elementLengths(match_direction) == 1) > 0){
					
					if("all" %in% directions | length(directions) > 1){
						print(paste("Found",sum(elementLengths(match_direction) == 1),"probes matching in",match_direction_name))
					}
					
					index <- which(elementLengths(match_direction) == 1)
					probeid <- rownames(probeData)[index]
					probeposition_here <- unlist(startIndex(match_direction)[index])
					names(probeposition_here) <- probeid
					positionVector <- c(positionVector,probeposition_here)
					
				}
			}
			return(positionVector)
		})		





setGeneric("plotOnGene", function(object,gene,probeData=NULL,label=NULL,genename=NULL,summaryType="median",interval=NULL,ylim=NULL,testType=NULL,forcePValue=FALSE,verbose=TRUE,cutoff=0.2,directions="all",ylab="expression") standardGeneric("plotOnGene"))
setMethod("plotOnGene", "ExpressionSet",
		function(object,gene,probeData=NULL,label=NULL,genename=NULL,summaryType="median",interval=NULL,ylim=NULL,testType=NULL,forcePValue=FALSE,verbose=TRUE,cutoff=0.2,directions="all",ylab="expression"){
			# This functions takes the following input:
			#  * "expressionset" the expression of all probes and their associated pData. Rownames of probeData and featurenames of expressionset should be the same.
			#  * "gene" in the form of a loaded FASTA sequence, character-sequence or DNAString instance
			#  * "probeData" (Optional if ProbeLevelSet is given) a description of all features. It requires a data frame with the featureNames as rownames and a column "sequence" containing their sequence.
			# and outputs a plot of the expression of individual probes as a function of location
			#
			# Optional switches are:
			# * label: instead plots the mean or median (see summarytype) of each of a set of factors, typically genotypes.
			# * summaryType: 
			#			options are:
			#				"mean" 
			#				"median"
			#				"quartiles"
			#				"dots", depending on which type of summary for each label is wanted.
			# * genename: Optional genename to put in title. Otherwise defaults to the title of the given FASTA sequence.
			# * interval: Optional interval of the gene to zoom on. If numbers outside of the gene length is given, they're truncated to end with the gene
			# * yMax: Optional max on Y-axis, used to zoom in.
			# * testType: Optional statistic procedure to identify especially interesting genes.
			#			options are:
			#				"linear model" - takes three factors (three genotypes for example) calculates the linear regression between them and circles significant value (at p<0.2) with increasing size circles depending on the significance.
			#				"students" - regular student's t-test between two factors
			#				"wilcoxons" - wilcoxon non-parametric test between two factors.
			# * forcePValue: Optional switch to force the writing of p-values
			# * verbose: Print a lot of outputtext?
			# * cutoff: at what p-value should probes be circled and printed. Defaults to 0.2. All probes with P >0.05 will be circled in grey, if cutoff is higher than 0.05.
			# * directions: A character vector of the matching-directions that should be scanned. Defaults to "all", which is shorthand for all the remaining directions: c("matchForwardSense","matchForwardAntisense","matchReverseSense","matchReverseAntisense")
			expressionset <- object
			#library(Biostrings)
			
			
			supported_summaryTypes <- c("mean","median","quartiles","dots")
			if(!summaryType %in% supported_summaryTypes)stop(paste("The summaryType",summaryType,"was not recognised. Only:",paste(supported_summaryTypes,collapse=", "),"is supported"))
			
			if(!(class(expressionset)[1] == "ExpressionSet" | class(expressionset)[1] == "ProbeLevelSet"))stop(paste("The class of the expressionset-variable, was not 'ExpressionSet' or 'ProbeLevelSet', but",class(expressionset)[1]," - this was not recognised with this function."))
			
			if(!is.null(label)){
				if(!(label %in% colnames(pData(expressionset)))){stop(paste("the label",label,"was not found in the pData of expressionset",expressionset@experimentData@title))}
				if(!(summaryType == "median"|summaryType == "mean"|summaryType == "quartiles"|summaryType == "dots")){stop(paste("the summaryType",summaryType,"is not recognised"))}
				label_class <- class(pData(expressionset)[, label])
				if(label_class %in% c("character","numeric","factor","integer")){
					if(label_class != "factor"){
						pData(expressionset)[, label] <- as.factor(pData(expressionset)[, label])
						if(verbose)print(paste("The label",label,"of class",label_class,"was temporarily converted to a factor with",length(levels(pData(expressionset)[, label])),"levels"))
					}
				}else{
					stop(paste("the label",label,"from",expressionset@experimentData@title,"was not of the type numeric, factor or character"))
				}
			}
			
			gene <- readGeneInput(gene,genename)
			if(length(gene) > 1){
				warning(paste("The input specified",length(gene),"genes, but the function only accepts one. The first one,",gene[[1]][["desc"]],"was used"))
				gene <- list(gene[[1]])
			}
			genename <- gene[[1]][["desc"]]
			sequence <- DNAString(gene[[1]][["seq"]])
			
			
			#This part choses if probe data is from the pData(featureData(x)) or from user-defined probe data (ie. a data frame with sequences)	
			if(is.null(probeData)){
				if(ncol(pData(featureData(expressionset))) == 0){
					stop(paste("No featuredata found in",expressionset@experimentData@title,"and no probe data explicitly given"))
				}else{
					probeData <- pData(featureData(expressionset))
					if(verbose)print(paste("Probe data for",nrow(probeData),"probes found as featureData in the expressionset",expressionset@experimentData@title))
				}
			}
			if(class(probeData[,"sequence"]) != "character")stop(paste("The 'sequence' column of the probeData had data of class",class(probeData[,"sequence"]),"and not character as required"))
			
			if(length(rownames(probeData)) == length(featureNames(expressionset))){
				if(!all(rownames(probeData) %in% featureNames(expressionset)))warning(paste("Not all probe ids given in the probe_level_annotation data in the notes section() of",expressionset@experimentData@title,"were found in the featureNames (probe ids) in the expressionset"))
				if(!all(featureNames(expressionset) %in% rownames(probeData)))warning(paste("Not all featureNames (probe ids) in the expressionset,",expressionset@experimentData@title,"were found in the probe ids given in the probe_level_annotation data in the notes section()"))
			}else{
				warning(paste("The expressionset,",expressionset@experimentData@title,"had",length(featureNames(expressionset)),"features in the exprs-part, but",length(rownames(probeData)),"entries in the probe_level_annotation of the notes() section. This might indicate corrupted data"))
			}
			
			
			#this part checks the interval
			if(!is.null(interval)){
				if(interval[1] > interval[2]){stop(paste("The interval must be a vector with two numeric values, the first being lowest. You gave:",paste(interval,collapse=" to ")))}
				if(interval[1]<1){
					interval[1] <- 1
					warning("Interval can't be lower than 1. It is automatically corrected to 1")
				}
				if(interval[2] > length(sequence)){
					warning(paste("Interval max (",interval[2],") can't be higher than gene length. It was automatically corrected to gene length ",length(gene),sep=""))
					interval[2] <- length(sequence)
				}
				sequence <- sequence[interval[1]:interval[2]]	
			}
			if(is.null(interval)){
				interval <- c(0,length(sequence))
			}
			interval_length <- interval[2]-interval[1]
			
			positionVector<-findProbePositions(object=expressionset, gene=list(list(desc=genename,seq=as.character(sequence))), probeData=probeData,interval=interval,directions="all",verbose=TRUE)
			
			if(length(positionVector) == 0){stop(paste("No probes found that matched in the given sequence of",genename,"- try specifying another direction variable with one of these: \"matchForwardSense\",\"matchForwardAntisense\",\"matchReverseSense\",\"matchReverseAntisense\""))}
			
			
			
			xlim <- c(interval[1],interval[2])
			if(is.null(ylim)){
				if(summaryType == "dots"){ylim <- c(0,max(exprs(expressionset[names(positionVector),])))	}
				if(summaryType == "mean"){ylim <- c(0,1.1*max(apply(exprs(expressionset[names(positionVector),]),1,mean)))}
				if(summaryType == "median"){ylim <- c(0,1.1*max(apply(exprs(expressionset[names(positionVector),]),1,median)))}
				if(summaryType == "quartiles"){ylim <- c(0,max(apply(exprs(expressionset[names(positionVector),]),1,quantile)["75%",]))}
				
#				ylim <- c(0,max(exprs(expressionset[names(positionVector)])))
			}
			
			
			
			#ylim <- c(0,1)
			##plot frame
			plot.default(1,
					type="n",
					xlab="bp",
					ylab=ylab,
					#yaxt="n",
					frame=FALSE,
					xlim=xlim,
					ylim=ylim,
					main=paste("Expression of probes"),
					cex.main=1
			)
			
			
			if(!is.null(label)){
				factors <- levels(pData(expressionset)[,label])
				oldc <- Sys.getlocale("LC_COLLATE")
				on.exit(Sys.setlocale("LC_COLLATE", oldc))
				Sys.setlocale("LC_COLLATE", "C")
				factors <- sort(factors)
				
			}
			if(is.null(label)){factors <- "one measurement"}
			
			if(length(factors) == 1){colours <- c("black")}
			if(length(factors) == 2){colours <- c("blue","red")}
			if(length(factors) > 2&length(factors)<5){colours <- c("blue","purple","red","yellow")}
			if(length(factors) >= 5&length(factors)<11){colours <- c("aquamarine","cyan","darkcyan","blue","blueviolet","darkred","red","orange","yellow","green")}
			if(length(factors) >= 11){
				stop("The use of labels with more than ten variables have not been implemented yet. Please work with the colours.")
			}
			
			
			
			#making the legend (third line): if there are less than 5 different factors we also calculate the amount of samples in each group
			if(!is.null(label)){		
				legend <- vector()
				for(i in 1:length(factors)){
					number_in_group <- sum(pData(expressionset)[,label] == factors[i],na.rm=TRUE)
					if(length(factors)<5){
						legend[i] <- paste(number_in_group,factors[i],"samples are",colours[i])
					}else{
						legend[i] <- paste(factors[i],"is",colours[i])
					}
				}	
				legend <- paste(legend,collapse=" - ")
			}
			
			
			mtext(genename,padj=-1.4,cex=0.7)
			
			if(!is.null(label)){
				if(summaryType == "median")mtext(paste("Each dot represents the median of probes from the following variables:",paste(factors,collapse=", "),"in",label),padj=-0.2,cex=0.7)
				if(summaryType == "mean")mtext(paste("Each dot represents the mean of probes from the variables:",paste(factors,collapse=", "),"in",label),padj=-0.2,cex=0.7)
				if(summaryType == "dots")mtext(paste("Each dot represents a sample. Colour codings from the following variables:",paste(factors,collapse=", "),"in",label),padj=-0.2,cex=0.7)
				if(summaryType == "quartiles")mtext(paste("Each line represents the 25% and 75% quartiles of probes from the following variables:",paste(factors,collapse=", "),"in",label),padj=-0.2,cex=0.7)
				mtext(legend,padj=1,cex=0.7)
			}
			
			if(is.null(label)){
				if(summaryType == "median")mtext("Each dot represents the median of probes at the given location",padj=-0.2,cex=0.7)
				if(summaryType == "mean")mtext("Each dot represents the mean of probes at the given location",padj=-0.2,cex=0.7)
				if(summaryType == "dots")mtext("Each dot represents a sample at the given location",padj=-0.2,cex=0.7)
				if(summaryType == "quartiles")mtext("Each line represents the 25% and 75% quartiles of probes at the given location",padj=-0.2,cex=0.7)
				#mtext(legend,padj=1,cex=0.7)
			}
			
			mtext(paste("The data is taken from the data set:",experimentData(expressionset)@title),padj=2.2,cex=0.7)
			
			
			y <- vector()
			for(i in 1:length(names(positionVector))){
				for(j in 1:length(factors)){
					x <- positionVector[i]+interval[1]
					
					if(!is.null(label)){expression_vector <- exprs(expressionset[names(positionVector)[i],])[pData(expressionset)[,label] == factors[j]]}
					if(is.null(label)){expression_vector <- exprs(expressionset[names(positionVector)[i],])}
					
					if(summaryType == "median"){y[j] <- median(expression_vector,na.rm=TRUE)}
					if(summaryType == "mean"){y[j] <- mean(expression_vector,na.rm=TRUE)}
					if(summaryType == "quartiles"){
						x <- positionVector[i]+interval[1]+(j*interval_length/480)-2   #the last terms makes sure that the lines are not on top of each other
						quantiles <- quantile(expression_vector,na.rm=TRUE)
						lines(x=c(x,x),y=c(quantiles[2],quantiles[3]),col=colours[j])
						lines(x=c(x-interval_length/300,x+interval_length/300),y=c(quantiles[3],quantiles[3]),col=colours[j])    #the ends of the line
						lines(x=c(x-interval_length/300,x+interval_length/300),y=c(quantiles[2],quantiles[2]),col=colours[j])    #the ends of the line
					}
					if(summaryType == "dots"){
						x <- positionVector[i]+interval[1]+(j*interval_length/480)-2   #the last terms makes sure that the dots are not on top of each other
						expression_vector <- expression_vector[is.finite(expression_vector)]
						points(x=rep(x,length(expression_vector)),y=expression_vector,col=colours[j],pch=19,cex=0.2)
					}
					if(summaryType == "mean"|summaryType == "median"){points(x=x,y=y[j],col=colours[j],pch=19,cex=1)} #plotting dots for mean and median summary types
				}
				if(summaryType == "mean"|summaryType == "median"){lines(x=c(x,x),y=c(min(y),max(y)))} #drawing lines for connecting mean and median dots
			}
			if(!is.null(testType)){
				plotStatistics(expressionset,probeData=probeData,label=label,summaryType=summaryType,testType=testType,interval=interval,forcePValue=forcePValue,verbose=verbose,positionVector=positionVector,ylim=ylim,xlim=xlim,cutoff=cutoff)
			}
		})



exonStructure <- function(mrna,genome,maxMismatch=4,y=0){
	# This program draws the exon structure on a gene
	# The genome variable is the sequence as downloaded through UCSC genome browser with exon structure splitted up.
	#
	# Importantly, the exon-numbers technically refer to "number of exon in investigated transcript", since fx in a download of exon 
	# structure for an isoform which skips exons, there will be skips of exon numbers.  
	
	mrna <- readGeneInput(mrna)
	genome <- readGeneInput(genome)
	
	exonstructure <- matrix(ncol=3,nrow=length(genome))
	exonrownames <- vector()
		
	if(length(mrna) > length(genome))stop(paste("There is",length(mrna),"entries in the mrna, and",length(genome),"entries in the genome argument. Are you sure you didn't switch the arguments?"))
	
	mrna <- DNAString(mrna[[1]][["seq"]][1])
	for(i in 1:length(genome)){
		exon <- DNAString(genome[[i]][["seq"]][1])
		match <- matchPattern(exon,mrna,max.mismatch=maxMismatch)
		if(length(start(match)) > 1){
			warning(paste("WARNING: exon",i,"was found to match",length(start(match)),"times on the mrna given. Only the first match have been used. This should be further investigated."))	
			exonstructure_here <- c(start(match)[1],width(match)[1]+start(match)[1],i)
		}else{
			exonstructure_here <- c(start(match),width(match)+start(match),i)
		}
		if(length(exonstructure_here) == 1){
			exonstructure_here <- c(NA,NA,NA)
			print(paste("Did not find a match in the mrna for exon",i,". If this is essential try using the exonStructure and increasing the maxMismatch number."))	
		}
		exonstructure[i,] <- exonstructure_here
	}

	colnames(exonstructure) <- c("start","end","exonnumber")
	
	for(i in 1:nrow(exonstructure)){
		if(!is.na(exonstructure[i,"start"])){
			lines(x=c(exonstructure[i,"start"],exonstructure[i,"end"]),y=c(y,y))
			points(x=exonstructure[i,"start"],y=y,pch="|")
			points(x=exonstructure[i,"end"],y=y,pch="|")
			text(x=mean(c(exonstructure[i,"start"],exonstructure[i,"end"])),y=y,labels=as.character(exonstructure[i,"exonnumber"]),cex=0.5,pos=1)
		}
	}
	
}


setGeneric("geneRegionScan", function(object,gene,genomicData=NULL,probeData=NULL,label=NULL,genename=NULL,summaryType="median",ylim=NULL,testType=NULL,forcePValue=FALSE,verbose=TRUE,cutoff=0.2,directions="all",correlationCutoff=0.3,probeLevelInfo=c("probeid")) standardGeneric("geneRegionScan"))
setMethod("geneRegionScan", "ExpressionSet",
		function(object,gene,genomicData=NULL,probeData=NULL,label=NULL,genename=NULL,summaryType="median",ylim=NULL,testType=NULL,forcePValue=FALSE,verbose=TRUE,cutoff=0.2,directions="all",correlationCutoff=0.3,probeLevelInfo=c("probeid")){
			#Wrapper around plotOnGene, exonStructure and plotCoexpression.
			#Serves a complete pdf-plot of the genes under investigation, their exon structure and their coexpression patterns.
			#The most important parameters are the "gene" which can be a vector of characters, DNAStrings or readFASTA-genes (it is parsed by readGeneInput) and the expressionset
			#Other tuning parameters are refered directly to the plotOnGene and plotCoexpression functions
			#See help file for these.
			#Output is named Report_ and the expressionsettitle and saved as pdf file in the working directory
			#Optional parameters, in addition to the ones for plotOnGene and plotCoexpression:
			#
			#	genomicData: This is  a list of readFASTA format of genomic data. It needs to be in the format downloaded from UCSC browser DNA-sequence downloader, using "one fasta record per region"
			#					and no introns. If specified, the plots will have exon-structure added to them.
			
			expressionset <- object
			
			gene <- readGeneInput(gene,genename)
			if(!is.null(genomicData)){
				if(length(gene) == 1){
					if(class(genomicData[[1]][[1]]) == "list"){ #a list(genomic) mistake
						genomicData <- genomicData[[1]]
					}
					genomicData <- readGeneInput(genomicData)
					genomicData <- list(genomicData)	
					
				}else{
					if(length(gene) != length(genomicData))stop("The genomicData and the gene did not have the same length. This is necessary")
					if(class(genomicData) != "list")stop(paste("The genomicData was not of class 'list', but of class",class(genomicData),"- this is not allowed when comparing multiple genes, and you are comparing",length(gene)))
					newGenomicData <- list()
					for(i in 1:length(genomicData)){
						newGenomicData[[i]] <- readGeneInput(genomicData[[i]])
					}
					genomicData <- newGenomicData
				}
			}
			
			if(verbose)print(paste("Now investigating",length(gene),"genes"))
			
			
			height=30
			
			pdf(width=length(gene)*(height/2),height=height,file=paste("report_",expressionset@experimentData@title,".pdf",sep=""))
			
			
			vertical_names <- split.screen(c(2,1))
			top_horizontal_names <- split.screen(c(1,length(gene)),screen=1)
			
			
			for(i in 1:length(gene)){
				screen(top_horizontal_names[i])
				gene_input <- list(gene[[i]])
				left_border  <-  1.43 * (length(gene) / 2)
				right_border  <-  1.02 * (length(gene) / 2)
				par(mai=c(0,left_border,0.82,right_border)) #this has been manually fine tuned to match the other functions called
				plotOnGene(expressionset,gene_input,probeData=probeData,label=label,genename=genename,summaryType=summaryType,ylim=ylim,testType=testType,forcePValue=forcePValue,verbose=verbose,cutoff=cutoff,directions=directions)
				if(!is.null(genomicData)){
					exonStructure(gene_input,genomicData[[i]])
				}
				
			}
			
			screen(vertical_names[2])
			plotCoexpression(expressionset,gene,probeData=probeData,verbose=verbose,directions=directions,correlationCutoff=correlationCutoff,probeLevelInfo=probeLevelInfo)
			close.screen(all.screens = TRUE)
			dev.off()
			
			
		})



readGeneInput <- function(gene,genename=NULL,verbose=TRUE){
#Function that will take a number of genes either as a vector of characters, a path to a fasta format file, as a vector of DNAstrings or as a readFASTA format. It will then output them in FASTA format for use with
#other programs. Optional argument genename forces a new name
	
	if(class(gene) == "character" & (length(grep("/",gene)) > 0 | length(grep("\\\\",gene,fixed=TRUE)) > 0)){
		if(file.exists(gene)){
			if(verbose)print(paste("Received what was interpreted as the path to a fasta format gene:",gene))
			#library(Biostrings)
			gene <- readFASTA_replacement(gene)
		}else{
			stop(paste("Received what was interpreted as the path to a fasta format gene, but did not find the file:",gene))
		}
	}
	
	if(class(gene) == "DNAString"){
		#library(Biostrings)
		gene <- as.character(gene)
	}
	
	if(class(gene) == "DNAStringSet"){
		gene <- as.character(gene)
	}
	
	
	if(class(gene) == "list"){
		if(class(gene[[1]]) == "DNAString"){#special case of a list of DNAStrings - probably only for deprecated usages of the Biostrings package
			new_gene <- vector()
			for(i in 1:length(gene)){
				new_gene <- c(new_gene,as.character(gene[[1]]))
			}
			gene <- new_gene
		}
	}
	
	
	if(class(gene) == "character"){
		new_gene <- list()
		for(i in 1:length(gene)){
			new_gene[[i]] <- list()
			if(!is.null(names(gene[i]))){
				new_gene[[i]][["desc"]] <- names(gene[i])
			}
			new_gene[[i]][["seq"]] <- gene[i]
		}
		
		gene <- new_gene
	}
	
	
	
	
	if(class(gene) == "list"){
		if(!"seq" %in% names(gene[[1]])){ #This would not be a fasta-format gene
			stop("The gene given was a list like a fasta-format gene is supposed to be after readFASTA(), but it lacks the \"seq\" part")			
		}
	}
	
	
	#adding desc-part
	if(!is.null(genename)){
		if(length(genename) != 1 & length(genename) != length(gene)){
			warning(paste("The length of the vector 'genename' is",length(genename),"and the number of genes given is",length(gene),"- I don't understand how to distribute the names to each gene, so they just got 'Unknown genename'"))
			genename <- "Unknown genename"
		}
		for(i in 1:length(gene)){
			if(length(gene) == length(genename))gene[[i]][["desc"]] <- genename[i]
			if(length(genename) == 1)gene[[i]][["desc"]] <- genename
		}
		
	}else{
		for(i in 1:length(gene)){
			if(is.null(gene[[i]][["desc"]])){
				gene[[i]][["desc"]] <- "Unknown genename"
				if(verbose)print("A gene was processed for which no name was found in fasta-format, and no name was explicitly given by 'genename' variable. The name 'Unknown genename' was assigned.")
			}
		}
	}
	
	if(class(gene) != "list")stop(paste("The format of the gene: \"",class(gene),"\" was not recognised at all.",sep=""))
	
	for(i in 1:length(gene)){
		if(!all(c("seq","desc") %in% names(gene[[i]]))){
			stop(paste("There was a problem with the gene number",i,"in the given gene. I'm afraid you'll have to open up the source if you want to understand, because this is really weird. 
									Alternatively, try to submit the sequence in a different format."))	
		}
	}
	
	return(gene)
}


setGeneric("plotCoexpression", function(object,gene,probeData=NULL,verbose=TRUE,directions="all",correlationCutoff=0.5,probeLevelInfo=c("probeid")) standardGeneric("plotCoexpression"))
setMethod("plotCoexpression", "ExpressionSet",
		
		function(object,gene,probeData=NULL,verbose=TRUE,directions="all",correlationCutoff=0.5,probeLevelInfo=c("probeid")){
			#Function to investigate co-expression of probes in a region. Takes a number of genes as DNAstring, vectors of DNAStrings, character-vectors or readFASTA outputs.
			#Takes an expressionset - which either has probe data (probe_name + sequence) specified in notes(expressionset)\$probe_level_annotation, or else and explicit data frame of this.
			#Optional parameters:
			# * directions: A character vector of the matching-directions that should be scanned. Defaults to matchForwardSense, 
			#		but can take fx: c("matchForwardSense","matchForwardAntisense","matchReverseSense","matchReverseAntisense")
			# * correlationCutoff: A number between 0 and 1, specifying at which cutoff correlations between probes should no longer be drawn. Defaults to 0.5
			# * probeLevelInfo: A character vector specifying what info should be supplied on plot for each probe. Options are probeid, probesetid, sequence. Default is only probeid.
			expressionset <- object
			
			debugging <- FALSE
			#library(Biostrings)
			
			
			
			supported_probeLevelInfo <- c("probeid","probesetid","sequence")
			if(!all(probeLevelInfo %in% supported_probeLevelInfo))stop(paste("Did not recognize all entries in probeLevelInfo. Choose only from:",paste(supported_probeLevelInfo,collapse=", ")))
			
			gene <- readGeneInput(gene)
			if(is.null(probeData)){
				if(ncol(pData(featureData(expressionset))) == 0){
					stop(paste("No featuredata found in",expressionset@experimentData@title,"and no probe data explicitly given"))
				}else{
					probeData <- pData(featureData(expressionset))
					if(verbose)print(paste("Probe data for",nrow(probeData),"probes found as featureData in the expressionset",expressionset@experimentData@title))
				}
			}
			if(class(probeData[,"sequence"]) != "character")stop(paste("The 'sequence' column of the probeData had data of class",class(probeData[,"sequence"]),"and not character as required"))
			
			
			if(length(rownames(probeData)) == length(featureNames(expressionset))){
				if(!all(rownames(probeData) %in% featureNames(expressionset))){
					warning(paste("Not all probe ids given in the probe_level_annotation data in the notes section() of",expressionset@experimentData@title,"were found in the featureNames (probe ids) in the expressionset"))
				}
				if(!all(featureNames(expressionset) %in% rownames(probeData))){
					warning(paste("Not all featureNames (probe ids) in the expressionset,",expressionset@experimentData@title,"were found in the probe ids given in the probe_level_annotation data in the notes section()"))
				}
			}else{
				warning(paste("The expressionset,",expressionset@experimentData@title,"had",length(featureNames(expressionset)),"features in the exprs-part, but",length(rownames(probeData)),"entries in the probe_level_annotation of the notes() section. This might indicate corrupted data"))
			}
			
#			probe <- XStringViews(probeData[,"sequence"],"DNAString")
			
			
			if(!all(directions %in% c("matchForwardSense","matchForwardAntisense","matchReverseSense","matchReverseAntisense","all"))){
				stop(paste("The directions variable:",paste(directions,collapse=", "),"was not recognised. The directions must be in \"matchForwardSense\",\"matchForwardAntisense\",\"matchReverseSense\",\"matchReverseAntisense\",\"all\""))
			}
			
			if("all" %in% directions){
				directions <- c("matchForwardSense","matchForwardAntisense","matchReverseSense","matchReverseAntisense")
			}
			
			
			positionVector <- vector()
			total_plot_size <- 1000
			
			padding_size <- 143 + (length(gene) * -8)   
			# the size of the gap between two gene sequences - will receive a vertical line if 
			#debugging =TRUE -- numbers were found by assuming linearity with length(gene) --- and 
			#then, yes, trial and error.
			
			if(debugging)print(paste("padding size is",padding_size))
			one_plot_size <- (total_plot_size - ((length(gene)-1)*( padding_size ))) / length(gene)
			if(debugging)x <- vector()
			for(i in 1:length(gene)){ 
				
				
				single_gene <- DNAString(gene[[i]][["seq"]])
				positionVector_here<-findProbePositions(object=expressionset, gene=single_gene, probeData=probeData,interval=interval,directions="all",verbose=TRUE)
				
#				if(length(positionVector) == 0){stop(paste("No probes found that matched in the given sequence of",genename,"- try specifying another direction variable with one of these: \"matchForwardSense\",\"matchForwardAntisense\",\"matchReverseSense\",\"matchReverseAntisense\""))}
				
				
#				positionVector_here <- vector()
#				probeid <- vector()
#				
#				
#				if("matchForwardSense" %in% directions){
#					probe_pdict <- PDict(probe)
#					matchForwardSense <- matchPDict(probe_pdict,single_gene)
#				}
#				if("matchForwardAntisense" %in% directions){
#					probe_pdict <- PDict(complement(probe))
#					matchForwardAntisense <- matchPDict(probe_pdict,single_gene)
#				}
#				if("matchReverseSense" %in% directions){ 
#					probe_pdict <- PDict(reverse(probe))
#					matchReverseSense <- matchPDict(probe_pdict,single_gene)
#				}
#				if("matchReverseAntisense" %in% directions){
#					probe_pdict <- PDict(reverseComplement(probe))
#					matchReverseAntisense <- matchPDict(probe_pdict,single_gene)
#				}
#				
#				for(match_direction_name in directions){
#					match_direction <- get(match_direction_name)
#					if(sum(elementLengths(match_direction) > 1) > 0){
#						stop(paste(sum(elementLengths(match_direction) > 1),"of the probes matched more than one time in the sequence"))
#					}
#					if(sum(elementLengths(match_direction) == 1) > 0){
#						index <- which(elementLengths(match_direction) == 1)
#						probeid <- rownames(probeData)[index]
#						probeposition_here <- unlist(startIndex(match_direction)[index])
#						names(probeposition_here) <- probeid
#						positionVector_here <- c(positionVector_here,probeposition_here)
#					}
#				}
				
				
				
				
				
				
				if(length(positionVector_here) == 0){print(paste("No probes found that matched in the given sequence"))}
				normalized_positionVector_here <- one_plot_size*(positionVector_here/length(single_gene))
				
				padding_size_here <- (i-1)*padding_size
				plots_before_size_here <- (i-1)*one_plot_size
				
				normalized_transposed_positionVector_here <- normalized_positionVector_here + plots_before_size_here + padding_size_here
				positionVector <- c(positionVector,normalized_transposed_positionVector_here)
				
				
				if(debugging){
					print(paste("The max position of gene",i,"was",round(max(positionVector_here)),"which was normalised and transposed to",round(max(normalized_transposed_positionVector_here)),"because the gene had length",length(single_gene)))
				}
				if(debugging){x <- c(x,plots_before_size_here+padding_size_here-padding_size,plots_before_size_here+padding_size_here)}
				if(debugging)print(paste("Will draw line at",plots_before_size_here+padding_size_here-padding_size,"and",plots_before_size_here+padding_size_here))
				
			}
			
			
			top_frame_size=total_plot_size * 0.04 * length(probeLevelInfo)
			if("sequence" %in% probeLevelInfo)top_frame_size <- top_frame_size*1.5
			frame()
			plot.window(xlim=c(0,total_plot_size),ylim=c(0,total_plot_size+top_frame_size))
			if(debugging)print(paste("making a window with size",paste(c(0,total_plot_size),collapse=", "),"on x, and",paste(c(0,total_plot_size+top_frame_size),collapse=", "),"on y"))
			
			if(debugging){
				x <- c(x,total_plot_size)
				for(x_here in x){
					
					lines(x=c(x_here,x_here),y=c(0,total_plot_size))
				}
			}
			
			
			
			if(debugging){
				print("class of of names(positionVector)") 
				print(class(names(positionVector)))
			}
			
			
			
			# This code block takes cares of the case where there are duplicated probe names (This should not really happen, but the program should be able to handle it even if it does)
			duplicates <- duplicated(names(positionVector))
			if(sum(duplicates) > 1 ){
				#this will need some extra processing - if this is not done the correlation between the same probe, but on different genes will really mess up. 
				if(verbose)print("Some probes matched twice in the given sequences. They have been marked with red colour")
				duplicated_probesets <- names(positionVector)[duplicates]
				duplicates <- names(positionVector) %in% names(positionVector)[duplicates]
				for(duplicated_probeset in duplicated_probesets){
					positions_of_duplicates <- names(positionVector) %in% duplicated_probeset
					for(i in 1:sum(positions_of_duplicates)){
						if(debugging)print(paste("renaming",names(positionVector[positions_of_duplicates][i]),"to",paste(names(positionVector[positions_of_duplicates][i]),i,sep="_")))
						names(positionVector)[positions_of_duplicates][i] <- paste(names(positionVector[positions_of_duplicates][i]),i,sep="_")
						if(i > 9){
							stop("There were more than 9 probeset duplicates (ie. probeset that matched more than one place in the given gene. This would cause the program to crash later. If necessary, it should be relatively easy to correct this later")
						}
					}
				}
				
			}
			
			
			#This code block writes the probeLevelInfo of the probes on screen
			size_to_write <- 0.5/ (length(probeLevelInfo) ^ 0.5)
			for(i in 1:length(positionVector)){
				text_to_write <- ""
				probeid <- names(positionVector[i])
				if(!probeid %in% featureNames(expressionset)){ #translating back to original probe name
					if(substr(probeid,nchar(probeid)-1,nchar(probeid)-1) == "_")probeid <- substr(probeid,0,nchar(probeid)-2)
				}
				
				if("probeid" %in% probeLevelInfo)text_to_write <- paste(text_to_write,probeid)
				if("probesetid" %in% probeLevelInfo)text_to_write <- paste(text_to_write,probeData[probeid,"probeset_name"])
				if("sequence" %in% probeLevelInfo)text_to_write <- paste(text_to_write,probeData[probeid,"sequence"])
				x = positionVector[i]
				if(duplicates[i]){
					colour <- "red"
				}else{
					colour <- "black"
				}
				text(x=x,y=total_plot_size*1.02,label=text_to_write,srt=-90,cex=size_to_write,adj=1,col=colour)
				lines(x=c(x,x),y=c(total_plot_size,total_plot_size*1.01))
			}
			
			
			
			
			#this code block calculates all the combinations and their correlation coefficents (also p-values for good measure) 
			combinations <- as.data.frame(t(combn(names(positionVector),2)))
			if(verbose& nrow(combinations) <= 2000)print(paste("Now calculating correlations for",nrow(combinations),"pairings of probes"))
			if(verbose & nrow(combinations) > 2000)print(paste("Now calculating correlations for",nrow(combinations),"pairings of probes - This might take a while..."))
			calculate_coexpression <- function(x){
				first_name <- x[1]
				second_name <- x[2]
				if(!first_name %in% featureNames(expressionset)|!second_name %in% featureNames(expressionset)){#translating back to original probe name
					if(substr(first_name,nchar(first_name)-1,nchar(first_name)-1) == "_")first_name <- substr(first_name,0,nchar(first_name)-2)
					if(substr(second_name,nchar(second_name)-1,nchar(second_name)-1) == "_")second_name <- substr(second_name,0,nchar(second_name)-2)
				}
				correlation <- cor.test(exprs(expressionset)[first_name,],exprs(expressionset)[second_name,])[["estimate"]]
				return(correlation)
			}
			
			correlation_calculations <- apply(combinations,1,calculate_coexpression)
			combinations <- cbind(combinations,correlation_calculations)
			colnames(combinations) <- c("first_name","second_name","correlation")
			combinations <- combinations[order(abs(combinations[,"correlation"])),]
			
			
			
			
			
			
			
			#making the scale
			positive_colour_scale <- brewer.pal(9,"Blues")
			negative_colour_scale <- brewer.pal(9,"Reds")
			negative_colour_scale <- negative_colour_scale[length(negative_colour_scale):1]
			colour_scale <- c(negative_colour_scale,positive_colour_scale)
			bins <- seq(-1,1,2/length(colour_scale))
			#scales are always the width of one_plot_size
			starting_point_of_scales <- length_of_scales_in_bits  <-  (total_plot_size *0.5) - (one_plot_size*0.5)
			for(i in 1:length(colour_scale)){
				x1 <- starting_point_of_scales + (i-1) * (one_plot_size / length(colour_scale))
				x2 <- starting_point_of_scales + (i) * (one_plot_size / length(colour_scale))		
				lines(x=c(x1,x2),y=c(0,0),col=colour_scale[i],lwd=4)
			}
			if(length(gene) > 1){
				text(x=starting_point_of_scales - one_plot_size * 0.1 ,y=total_plot_size/25,label="Correlation coefficent:",cex=0.7,adj=1)
				text(x=starting_point_of_scales - one_plot_size * 0.1 ,y=total_plot_size/50,label=paste("Cutoff: abs(coefficent) >",correlationCutoff),cex=0.7,adj=1)
			}else{
				text(x=starting_point_of_scales + one_plot_size * 0.1 ,y=total_plot_size/12,label="Correlation coefficent:",cex=0.7,adj=0)
				text(x=starting_point_of_scales + one_plot_size * 0.1 ,y=total_plot_size/25,label=paste("Cutoff: abs(coefficent) >",correlationCutoff),cex=0.7,adj=0)
				
			}
			for(i in seq(-1,1,0.5)){
				x <- total_plot_size * 0.5  + i * one_plot_size * 0.5
				text(x=x,y=total_plot_size/50,label=i,cex=0.7)	
			}
			
			
			
			for(i in 1:(nrow(combinations))){
				first_name <- as.character(combinations[i,"first_name"])
				second_name <- as.character(combinations[i,"second_name"])
				first_pos <- positionVector[first_name]
				second_pos <- positionVector[second_name]
				#print(paste(first_name,"at",round(first_pos),"and",second_name,"at",round(second_pos),"had correlation",combinations[i,"correlation"]))
				bin_no_here <- sum(bins<combinations[i,"correlation"])
				
				if(abs(combinations[i,"correlation"]) > correlationCutoff){ 
					#print(paste(first_name,"and",second_name,"with correlation",combinations[i,"correlation"],"was put in bin no",bin_no_here))
					lines(x=c(first_pos,mean(c(first_pos,second_pos))),y=c(total_plot_size,total_plot_size-abs(first_pos-second_pos)),col=colour_scale[bin_no_here],lty=2)
					lines(x=c(second_pos,mean(c(first_pos,second_pos))),y=c(total_plot_size,total_plot_size-abs(first_pos-second_pos)),col=colour_scale[bin_no_here],lty=2)
					points(x=mean(c(first_pos,second_pos)),y=total_plot_size-abs(first_pos-second_pos),col=colour_scale[bin_no_here],cex=abs(combinations[i,"correlation"])*2,pch=18)
					
					if(debugging){
						if(abs(first_pos-second_pos) == 0){
							print(paste(first_name,"at",round(first_pos),"of class",class(first_name),"and",second_name,"at",round(second_pos),"of class",class(second_name),"had correlation",combinations[i,"correlation"],"- this is weird because it is very high on the y-axis"))
							print(paste("The sequences are:",probeData[as.character(first_name),"sequence"],"and",probeData[as.character(second_name),"sequence"],"respectively"))
						}
					}
				}
			}
		})



