

setGeneric("translateSampleNames", function(object, translationFile, from, to) standardGeneric("translateSampleNames"))
setMethod("translateSampleNames", "ExpressionSet",
		function(object, translationFile, from, to) {
			# Function that can translate the sampleNames of an expressionset given a translation file.
			# The translationFile gives the translation between samplenames now and the samplenames as they should be. 
			# It is either the path to a tab-seperated txt-file with headers or a directly a data frame
			# The variables from and to specify the column names to translate.
			# They must be present in the translationFile and the from must be present in the sampleNames of the expressionset
			# Extra entries in the translationFile are omitted but samples in the expressionset without samples in the translationFile
			# will raise an error.
			expressionset <- object
			
			
			if(length(translationFile) == 1){
				if(file.exists(translationFile)){
					translationFile <- read.table(translationFile, sep = "\t", header = TRUE)
				}
			}
			if(class(translationFile) != "data.frame")stop(paste("The class of the translationFile -",class(translationFile),"- was not recognised either as a data.frame or a path to a tab-separated text file"))
			
			possibilites<-paste(colnames(translationFile),collapse=", ")
			if(missing(from)){
				stop(paste("The 'from' argument is required. It should specify one of the following columns in the translationFile:",possibilites))
			}
			
			
			if(! from %in% colnames(translationFile)){
				stop(paste("The translation file should contain the columnname",from,"as specified in the variable 'from'. The translation file only has the following:",possibilites))
			}
			
			if(missing(to)){
				stop(paste("The 'to' argument is required. It should specify one of the following columns in the translationFile:",possibilites))
			}
			
			if(! to %in% colnames(translationFile)){
				stop(paste("The translation file should contain the columnname",to,"as specified in the variable 'to'. The translation file only has the following:",possibilites))
			}
			
			if(! all(translationFile[,from] %in% sampleNames(expressionset))){
				
				missing<-translationFile[,from][! translationFile[,from] %in% sampleNames(expressionset)]
				
				stop(paste(length(missing),"samples in the given expressionset are not present in the",from,"column of the translationFile:",paste(missing,collapse=", ")))
			}
			
			
			expressionset_from_IDs <- sampleNames(expressionset)
			expressionset_to_IDs <- vector()
			
			for(expressionset_CEL_ID in expressionset_from_IDs){
				sorting_vector <- translationFile[,from]%in%expressionset_CEL_ID
				if(sum(sorting_vector) != 1)stop(paste("Major and really weird error with",expressionset_CEL_ID))
				expressionset_to_ID <- as.character(translationFile[sorting_vector,to])
				
				expressionset_to_IDs <- c(expressionset_to_IDs,expressionset_to_ID)
			}
			if(length(expressionset_to_IDs)  !=  length(unique(expressionset_to_IDs))){
				stop(paste("There is",length(expressionset_to_IDs) - length(unique(expressionset_to_IDs)),"non-unique new sample IDs after processing the expressionset as specified. This is not allowed"))
			}
			
			rownames(pData(expressionset)) <- expressionset_to_IDs
			colnames(exprs(expressionset)) <- expressionset_to_IDs
			
			return(expressionset)
		})



setGeneric("excludeDoubleMatchingProbes", function(object,genome="BSgenome.Hsapiens.UCSC.hg18",verbose=TRUE,directions=c("matchForwardSense","matchForwardAntisense","matchReverseSense","matchReverseAntisense"),previousData = NULL) standardGeneric("excludeDoubleMatchingProbes"))
setMethod("excludeDoubleMatchingProbes", "ProbeLevelSet",
		function(object,genome="BSgenome.Hsapiens.UCSC.hg18",verbose=TRUE,directions=c("matchForwardSense","matchForwardAntisense","matchReverseSense","matchReverseAntisense"),previousData = NULL){
			# Function that will automatically remove those probes in a ProbeLevelSet that have matches in more than one place in the genome
			
			
			if(is.null(previousData)){
				matching_result <- findSequenceInGenome(sequences=getSequence(object),genome=genome,verbose=verbose,directions=directions)
			}else{
				if(verbose)print("Received previously scanned data, and will use this instead")
				if(any(!as.character(previousData[,"sequence_name"]) %in% featureNames(object))){
					badprobes <- unique(as.character(previousData[,"sequence_name"])[!as.character(previousData[,"sequence_name"]) %in% featureNames(object)])
					if(length(badprobes)<10){
						stop(paste("The following probes in the given previousData was not found in the object:",paste(badprobes,collapse=", ")))
					}else{
						stop(paste(length(badprobes),"probes in the given previousData was not found in the object"))
					}
				}
				matching_result <- previousData
				
			}
			
			
			duplicated_logical <- duplicated(matching_result[,"sequence"])
			duplicated_sequences <- unique(matching_result[duplicated_logical,"sequence"])
			duplicated_probeids <- names(getSequence(object))[getSequence(object) %in% duplicated_sequences]
			if(length(duplicated_probeids) > 0){
				if(verbose)print(paste("Found",length(duplicated_probeids),"- they have been removed and output from genome scan have been saved in notes of expressionset"))
				keep_in_feature_data <- !rownames(pData(featureData(object))) %in% duplicated_probeids
				keep_in_exprs <- !rownames(exprs(object)) %in% duplicated_probeids
				
				pData(featureData(object)) <- pData(featureData(object))[keep_in_feature_data,]
				exprs(object) <- exprs(object)[keep_in_exprs,]
				
			}else{
				if(verbose)print(paste("Did not find any duplicated probes. Output from genome scan have been saved in notes of expressionset"))
			}
			notes(object)[["duplicates_in_genome_scan"]] <- matching_result
			return(object)
		})


findSequenceInGenome <- function(sequences,genome="BSgenome.Hsapiens.UCSC.hg19",verbose=TRUE,directions=c("matchForwardSense","matchForwardAntisense","matchReverseSense","matchReverseAntisense")){
	#This functions takes a set of sequences and checks where they are present in UCSC genome in the relevant bioconductor package.
	#The sequence should be given as a list of nucleotide-character strings: c("GAGTATA","GTAGATGA")
	#Optional arguments:
	#	genome: the name of the BSgenome package. Defaults to the human genome.
	#	verbose: controlling the amount of output
	#	directions: which directions (reverse, complement, reverse-complement, forward) to scan in defaults to all four: ("matchForwardSense",
	#										"matchForwardAntisense","matchReverseSense","matchReverseAntisense")
	#
	# This function will take quite a lot of time - depending on number of sequences to scan. 
	if(!require(BSgenome)) stop("Package BSgenome is required and not installed.")
	
	if(!all(directions %in% c("matchForwardSense","matchForwardAntisense","matchReverseSense","matchReverseAntisense"))){
		stop(paste("The directions variable:",paste(directions,collapse=", "),"was not recognised. The directions must be in \"matchForwardSense\",\"matchForwardAntisense\",\"matchReverseSense\",\"matchReverseAntisense\""))
	}
	if(verbose)print(paste("Now scanning the genome in the following directions:",paste(directions,collapse=", ")))
	if(!require(genome,character.only=TRUE))stop(paste("Did not find the package for the genome:",genome))
	
	
	genome_name <- sub("BSgenome.","",genome)
	genome_name <- sub("\\..+$","",genome_name)
	chr_names <- seqnames(get(genome_name))
	chr_names <- chr_names[-grep("_", chr_names)]
	
	
	
#	chr_names <- chr_names[18] #for testing
	
	if(class(sequences) == "list"){
		print("The sequences given was in the form of a list. These have been converted to a character vector")
		new_sequences <- unlist(lapply(sequences, as.character))
	}
	
	probe <- DNAStringSet(as.character(sequences)) #need to remove names or PDict will crash
	
	probe_directions <- DNAStringSet()
	if("matchForwardSense" %in% directions){
		probe_directions <- c(probe_directions,probe)
	}
	if("matchForwardAntisense" %in% directions){
		probe_directions <- c(probe_directions,complement(probe))
	}
	if("matchReverseSense" %in% directions){ 
		probe_directions <- c(probe_directions,reverse(probe))
	}
	if("matchReverseAntisense" %in% directions){
		probe_directions <- c(probe_directions,reverseComplement(probe))			
	}
	
	probe_pdict <- PDict(probe_directions)
	
	
	result <- data.frame()
	
	for(chr_name in chr_names){
		if(verbose){print(paste("scanning",chr_name,"for matches"))}
		
		chr_seq <- unmasked(get(genome_name)[[chr_name]])		
		
		mindex <- matchPDict(probe_pdict,chr_seq)
		positions <- startIndex(mindex)
		
		match_indices<-which(!unlist(lapply(positions,is.null)))
		
		entrynumber <- vector()
		hitposInChr <- vector()
		chr <- vector()
		sequence <- vector()
		
		for(match_index in match_indices){
			for(i in 1:length(positions[[match_index]])){
				# correcting for extended PDict to get correct the entrynumber (if reverse or complement was also included)				
				entrynumber_here<-match_index%%length(sequences)
				if(entrynumber_here == 0) entrynumber_here <- length(sequences)
				
				entrynumber <- c(entrynumber, entrynumber_here)
				hitposInChr  <- c(hitposInChr, positions[[match_index]][i])
				chr <- c(chr,chr_name)
				sequence <- c(sequence, sequences[entrynumber_here])
			}
		}
		
		
		if(verbose){print(paste("found",length(hitposInChr),"matches between a probe and",chr_name))}		
		result_this_chromosome <- data.frame(entrynumber=entrynumber,hitposInChr=hitposInChr,chr=chr,sequence=sequence)
		
		result <- rbind(result,result_this_chromosome)		
		rm(chr_seq)
		gc()
	}
	
	
	
	
	
	
	if(nrow(result) > 0){
		
		colnames(result) <- c("entrynumber","hitposInChr","chr","sequence")
		
		if(!is.null(names(sequences))){
			if(verbose)print("Found names for each of the sequences. These will be included in the output")
			sequence_name <- vector()
			
			
			for(number in result[,"entrynumber"]){
				sequence_name <- c(sequence_name,names(sequences)[number])
			}
			
			result <- cbind(result,sequence_name)
		}
		return(result)	
	}
	else{
		print("No matches were found in the genome")
		return(NULL)
	}
}



setGeneric("addSnpPdata", function(object,listOfSnps,individualNames="sampleNames") standardGeneric("addSnpPdata"))
setMethod("addSnpPdata", "ExpressionSet",
		function(object,listOfSnps,individualNames="sampleNames"){
			# takes an expressionset and the path for a SNP-file. The SNP-file should have the same format as 
			# files downloaded from hapmap.org
			# The function then checks if all samplenames are present in the SNP-file and vice versa. If they are it
			# decorates the expressionset with the new SNPs. It also checks for existing entries with the same name
			# and if they are found it adds the extra SNPs to the this entry.
			# optional argument individualNames can refer to a column in the pData to get individualIds (=whatever the listOfSnps call the individuals). Defaults to just sampleNames
			
			expressionset <- object
			
			if(individualNames == "sampleNames"){
				individualNamesVector <- sampleNames(expressionset)
			}else{
				if(!individualNames %in% colnames(pData(expressionset))) stop(paste("The individualNames",individualNames,"was not found as a column name in the pData of the expressionset"))
				individualNamesVector <- pData(expressionset)[,individualNames]
			}
			names(individualNamesVector) <- sampleNames(expressionset)
			
			pdata_raw <- read.table(listOfSnps,skip=2, comment.char = "",header=TRUE,row.names=1)
			snp_data <- pdata_raw[,c("alleles","chrom","pos","strand","assembly.","center","protLSID","assayLSID","panelLSID","QCcode")]
			pdata_raw <- t(pdata_raw[,!colnames(pdata_raw) %in% colnames(snp_data)])
			
			
			
			
			if(!any(rownames(pdata_raw) %in% individualNamesVector))stop("None of the sample names in the listOfSnps file was found in the expressionset")
			print(paste(sum(rownames(pdata_raw) %in% individualNamesVector),"of the",nrow(pdata_raw),"samplenames in the listOfSnps was found in the expressionset"))	
			
			
			for(snp in colnames(pdata_raw)){
				genotype_vector <- vector()
				if(snp %in% colnames(pData(expressionset))){
					genotype_first <- expressionset[[snp]]
					for(sampleName in names(individualNamesVector)){
						individualName <- individualNamesVector[sampleName]
						if(individualName %in% rownames(pdata_raw)){
							if(!as.character(pdata_raw[individualName,snp]) == "NN"){
								if(is.na(pData(expressionset)[sampleName,snp])){
									genotype_vector <- c(genotype_vector,as.character(pdata_raw[individualName,snp]))	
								}else{
									if(pData(expressionset)[sampleName,snp]  !=  pdata_raw[individualName,snp]){
										entry_1<-as.character(pData(expressionset)[sampleName,snp])
										entry_2<-as.character(pdata_raw[individualName,snp])
										print(paste("WARNING: the expressionset already contained the entry",entry_1,"for",snp,"and",individualName,"but the listOfSnps had",entry_2,"- the original has been preserved"))
									}
									genotype_vector <- c(genotype_vector,as.character(pData(expressionset)[sampleName,snp]))
								}
							}else{
								genotype_vector <- c(genotype_vector,as.character(pData(expressionset)[sampleName,snp]))
							}
						}else{
							genotype_vector <- c(genotype_vector,as.character(pData(expressionset)[sampleName,snp]))
						}
					}
					expressionset[[snp]] <- as.factor(genotype_vector)
					genotype_last <- expressionset[[snp]]
					
					if(sum(is.na(genotype_last)) != sum(is.na(genotype_first))){
						first_num<-sum(is.na(genotype_first))
						last_num<-sum(is.na(genotype_last))
						print(paste(snp,"was already found in the given expressionset. It had",first_num,"NA entries to begin with. Now it has",last_num))
					}
				}else{
					for(sampleName in names(individualNamesVector)){
						individualName <- as.character(individualNamesVector[sampleName])
						if(individualName %in% rownames(pdata_raw)){
							if(!as.character(pdata_raw[individualName,snp]) == "NN"){
								genotype_vector <- c(genotype_vector,as.character(pdata_raw[individualName,snp]))
							}else{
								genotype_vector <- c(genotype_vector,NA)
							}
						}else{
							genotype_vector <- c(genotype_vector,NA)
						}
					}
					pData(expressionset) <- cbind(pData(expressionset),as.factor(genotype_vector))
					colnames(pData(expressionset))[ncol(pData(expressionset))] <- snp
				}
			}
			
			notes(expressionset)[["snp_data"]] <- snp_data
			
			
			
			return(expressionset)
		})



getServerProbeIntensities <- function(
		listOfProbesets,
		celfilePath,
		annotation=NULL,
		aptCelExtractPath=NULL,
		pgfPath=NULL,
		clfPath=NULL,
		cdfPath=NULL,
		serveradress="localhost",
		username=NULL,
		password=NULL,
		plinkPath=NULL,
		pscpPath=NULL,
		verbose=TRUE){
#Wrapper around getLocalProbeIntensities, that will send all necesarry stuff to a server and perform the
# getLocalProbeIntensities calculations there
#
#	serveradress: the ip-adress of a remote server on which the script can be run
#	username: needed with serveradress specified
#	password: needed with serveradress specified
#
#
	local_working_dir <- getwd()
	remote_working_dir <- celfilePath
	
	#some checking of the input
	if(substr(celfilePath,1,1) == "~"){
		stop("Do not use ~ as short for your home directory, as this will cause problems in the interface between different operating systems in this function")
	}
	endChar <- substr(celfilePath,nchar(celfilePath),nchar(celfilePath))
	if(!endChar %in% c("/","\\")){
		paste("Because this function must work between different operating systems, you will have to make sure that the celfilePath given ends with the file separator of the remote system (ie. '/' or '\\')")
	}
	
	
	#The following block is because I want users to both be able to submit vectors of probesets, and paths for files with probesets
	if(class(listOfProbesets[1]) == "numeric"|class(listOfProbesets[1]) == "integer")listOfProbesets <- as.character(listOfProbesets)
	if(length(listOfProbesets) == 1 & file.exists(listOfProbesets[1])){
		delete_probes_file <- FALSE
		if(!(readLines(listOfProbesets,n=2,warn=FALSE)[1] %in% c("probeset_id","probeset_name"))){
			oldlines <- readLines(listOfProbesets)
			probe_list_file <- file(listOfProbesets,"w")
			writeLines(c("probeset_id",oldlines),probe_list_file)
			close(probe_list_file)
		}
	}else{
		if(length(listOfProbesets) == 1){
			print("WARNING: the list of probesets only had one entry. This entry will be treated as a name of a probeset, but it is also possible that it is a missing file")
		}
		listOfProbesets <- c("probeset_id",listOfProbesets)
		probe_list_file <- file("probelist","w")
		writeLines(listOfProbesets,probe_list_file)
		close(probe_list_file)
		listOfProbesets <- paste(getwd(),.Platform[["file.sep"]],"probelist",sep="")
		delete_probes_file <- TRUE
	}
	
	
	if(serveradress == "localhost"){
		#easy wrap
		expressionset <- getLocalProbeIntensities(
				listOfProbesets=listOfProbesets,
				celfilePath=celfilePath,
				annotation=annotation,
				aptCelExtractPath=aptCelExtractPath,
				pgfPath=pgfPath,
				clfPath=clfPath,
				cdfPath=cdfPath,
				verbose=verbose)
		
	}else{
		#Running on server check username, password, plink, pscp
		
		#if both plinkPath and pscpPath is set to NULL, we'll try to see if they are stored in the locations.txt file.
		locationFilePath <- file.path(.path.package("GeneRegionScan"),"configfiles","locations.txt")
		if(is.null(plinkPath) & is.null(pscpPath) & (file.access(locationFilePath,4) == 0)){
			locationFile <- file(locationFilePath,"r")
			locationFileLines <- readLines(locationFile)
			close(locationFile)
			splitLocationFileLines <- strsplit(locationFileLines,"\t")
			for(line in splitLocationFileLines){
				if("plinkPath" == line[1])plinkPath <- line[3]
				if("pscpPath" == line[1])pscpPath <- line[3]
			}
			#if(!is.null(plinkPath) & verbose)print(paste("plinkPath given as null, but plink has been used before and was set to:",plinkPath))
			#if(!is.null(pscpPath) & verbose)print(paste("pscpPath given as null, but pscp has been used before and was set to:",pscpPath))
		}
		
		if(is.null(plinkPath))stop("When specifying a serveradress, a plinkPath is also needed")
		if(is.null(pscpPath))stop("When specifying a serveradress, a pscpPath is also needed")
		if(is.null(username) | is.null(password)){stop("When running on remote server, you need to specify the variables username and password")}
		if(nchar(password)<1)stop("You must supply a password for putty to work - \"\" does not work")
		if(!file.exists(plinkPath) | !file.exists(pscpPath)){stop("One or more of the required programs plink or pscp was not found")}
		#if(!exists("getLocalProbeIntensities")){stop("getLocalProbeIntensities was not found. This must be available in the environment")}
		#if(!exists("getProbeLevelAnnotationForExonArrays")){stop("getProbeLevelAnnotationForExonArrays was not found. This must be available in the environment")}
		
		transfer_file <- function(remote_path,local_path,type){
			if(type == "retrieve")system(paste(shQuote(pscpPath)," -pw ",password," ",username,"@",serveradress,":",shQuote(remote_path)," ",shQuote(local_path),sep=""))
			if(type == "send")system(paste(shQuote(pscpPath)," -pw ",password," ",shQuote(local_path)," ",username,"@",serveradress,":",shQuote(remote_path),sep=""))
		}
		
		
		save(list = c("annotation","celfilePath","aptCelExtractPath","clfPath","pgfPath","cdfPath"), file = "shipment_to_server.rdata")
		transfer_package_name_and_path <- file.path(getwd(),"shipment_to_server.rdata")
		transfer_file(paste(celfilePath,"shipment_to_server.rdata",sep=""),transfer_package_name_and_path,"send")
		
		
		transfer_file(paste(celfilePath,"probelist",sep=""),listOfProbesets,"send")
		
		
		
		server_command_file_path_and_name <- file.path(getwd(),"server_file.txt")
		server_command_file <- file(server_command_file_path_and_name,"w")
		writeLines(c(paste("R -f ",celfilePath,"r_file.txt --no-save",sep="")),server_command_file)
		close(server_command_file)
		
		r_command_file_path_and_name <- file.path(getwd(),"r_file.txt")
		r_command_file <- file(r_command_file_path_and_name,"w")
		writeLines(c(
						paste("setwd(\"",celfilePath,"\")",sep=""),
						"if(!require(GeneRegionScan))stop(\"GeneRegionScan package must also be installed on remote computer!\")",
						"load(\"shipment_to_server.rdata\")",
						"expressionset<-getLocalProbeIntensities(paste(celfilePath,\"probelist\",sep=\"\"),celfilePath,annotation,aptCelExtractPath,pgfPath,clfPath,cdfPath)",
						"rm(list=ls()[ls() != \"expressionset\"])",
						"save.image(\"newexpressionset.rdata\")",
						"quit(save=\"no\")"
				),r_command_file)
		close(r_command_file)
		transfer_file(paste(celfilePath,"r_file.txt",sep=""),r_command_file_path_and_name,"send")
		
		system(paste(shQuote(plinkPath)," ",username,"@",serveradress," -pw ",password," -m ",shQuote(server_command_file_path_and_name),sep=""))
		
		transfer_file(paste(celfilePath,"newexpressionset.rdata",sep=""),"newexpressionset.rdata","retrieve")
		
		server_clean_file_path_and_name <- file.path(getwd(),"server_file.txt")
		server_clean_file <- file(server_clean_file_path_and_name,"w")
		writeLines(c(
						paste("rm ",celfilePath,"newexpressionset.rdata",sep=""),
						paste("rm ",celfilePath,"cellist",sep=""),
						paste("rm ",celfilePath,"probelist",sep=""),
						paste("rm ",celfilePath,"shipment_to_server.rdata",sep=""),
						paste("rm ",celfilePath,"r_file.txt",sep="")),
				server_clean_file)
		
		close(server_clean_file)
		system(paste(shQuote(plinkPath)," ",username,"@",serveradress," -pw ",password," -m ",shQuote(server_clean_file_path_and_name),sep=""))
		
		unlink(server_clean_file_path_and_name)
		unlink(r_command_file_path_and_name)
		unlink(server_command_file_path_and_name)
		unlink("shipment_to_server.rdata")
		load("newexpressionset.rdata")
		if(!exists("expressionset"))stop("The completed expressionset was not found. This is a serious error")
		unlink("newexpressionset.rdata")
		
		#saving location of plink pscp for the next time
		#this functionality has been removed because it is not allowed to save settings between sessions in Bioconductor.
		#it was re-inserted because other it seems that other packages do it as well 
		if(file.access(locationFilePath,2) == 0){
			locationFile <- file(locationFilePath,"r")
			locationFileLines <- readLines(locationFile)
			close(locationFile)
			splitLocationFileLines <- strsplit(locationFileLines,"\t")
			newLocationFileLines <- vector()
			for(line in splitLocationFileLines){
				if(!line[1] %in% c("plinkPath","pscpPath")){
					newLocationFileLines <- c(newLocationFileLines,paste(line[1],line[2],line[3],sep="\t"))
				}
			}	
			newLocationFileLines <- c(newLocationFileLines,paste("plinkPath","all",plinkPath,sep="\t"))
			newLocationFileLines <- c(newLocationFileLines,paste("pscpPath","all",pscpPath,sep="\t"))
			locationFile <- file(locationFilePath,"w")
			writeLines(newLocationFileLines,con=locationFile)
			close(locationFile)
			
		}
		if(file.access(locationFilePath,2) != 0 & verbose)print(paste("Would have saved the locations of plinkPath and pscpPath for next time, but did not have write 
									permission for the configfile",locationFilePath))
		
		
		
		
	}
	if(delete_probes_file)unlink(listOfProbesets)
	return(expressionset)
}


getLocalProbeIntensities <- function(listOfProbesets,celfilePath,annotation=NULL,aptCelExtractPath=NULL,pgfPath=NULL,clfPath=NULL,cdfPath=NULL,verbose=TRUE){
#function that will take a list of probe sets and a celfilePath. And retrive the expression level of the probes in the
#probe set. These will be made into an expressionset, in which the sequence of each probe is also included in the notes.
#It is suggested that this set is run through the verify_position afterwards to expand the info given.
#
#Arguments are:
#	listOfProbesets:	Either a vector of probe names or the path and name of a file containing the probe names
#	celfilePath: pathname to a folder with cel files for the given dataset
#	annotation: string with the annotation. Can be HuEx-1_0-st-v2 or hgu133plus2 at the moment.
#	
	if(is.null(cdfPath) & is.null(pgfPath) & is.null(clfPath) & is.null(annotation)){
		stop("You need to provide at least one of the following: pgfPath, clfPath, cdfPath or annotation - and probably more, but error messages will tell you which")
	}
	if(!is.null(cdfPath) & (!is.null(pgfPath) | !is.null(clfPath))){
		stop("Do not specify both cdfPath and a clfPath or pgfPath. pgfPath and clfPath is for exon arrays that do not use cdf files")
	}
	
	if(!file.exists(celfilePath))stop(paste("The celfilePath",celfilePath,"was not found"))
	if(!file.info(celfilePath)[["isdir"]])stop(paste("The celfilePath",celfilePath,"was not identified as a directory. You must specify a directory."))
#	if(substr(celfilePath,nchar(celfilePath),nchar(celfilePath)) == .Platform[["file.sep"]]){
#		celfilePath <- substr(celfilePath,1,nchar(celfilePath)-length(.Platform[["file.sep"]]))
#	}
	celfilePath<-paste(dirname(celfilePath),basename(celfilePath),sep=.Platform[["file.sep"]])
	
	celfilePath <- path.expand(celfilePath)
	
	if(!is.null(annotation)){
		if(length(grep("pgfPath",annotation)) + length(grep("clfPath",annotation)) + length(grep("cdfPath",annotation)) > 0){
			stop("For some obscure reason the words cdfPath, pgfPath or clfPath appears in the annotation name. This will cause the program to crash later. Please change the annotation name.")
		}
	}
	
	#if all cdf, pgf and clf is set to null, but annotation is not, we'll assume it is a lazy user
	#and see if this type of annotation has been used before.
	locationFilePath <- file.path(.path.package("GeneRegionScan"),"configfiles","locations.txt")
	if(is.null(cdfPath) & is.null(pgfPath) & is.null(clfPath) & !is.null(annotation) & (file.access(locationFilePath,4) == 0)){
		locationFile <- file(locationFilePath,"r")
		locationFileLines <- readLines(locationFile)
		close(locationFile)
		splitLocationFileLines <- strsplit(locationFileLines,"\t")
		for(line in splitLocationFileLines){
			if(annotation == line[2]){
				if("pgfPath" == line[1])pgfPath <- line[3]
				if("clfPath" == line[1])clfPath <- line[3]
				if("cdfPath" == line[1])cdfPath <- line[3]
			}
		}
		if(!is.null(cdfPath) & verbose)print(paste("cdfPath given as null, but the annotation",annotation,"has used this file before:",cdfPath))
		if(!is.null(clfPath) & verbose)print(paste("clfPath given as null, but the annotation",annotation,"has used this file before:",clfPath))
		if(!is.null(pgfPath) & verbose)print(paste("pgfPath given as null, but the annotation",annotation,"has used this file before:",pgfPath))
	}
	#
	annotation_type <- ""
	if(is.null(cdfPath) & is.null(pgfPath) & is.null(clfPath))stop(paste("For the annotation",annotation,"you need to provide at least one of the following: pgfPath, clfPath, cdfPath"))
	if(is.null(cdfPath)){
		#assuming exon type array
		if(is.null(pgfPath)|is.null(clfPath))stop("This annotation needs both a pgfPath argument and a cdfPath argument")
		if(!file.exists(pgfPath)|!file.exists(clfPath))stop(paste("The pgf or clf file was not found"))
		if(file.info(pgfPath)[["isdir"]]|file.info(clfPath)[["isdir"]])stop(paste("The pgf or clf specified is a directory, not a file"))
		annotation_type <- "doesnothavecdf"
	}
	if(is.null(clfPath)|is.null(pgfPath)){
		#assuming 3'IVT type array (ie with already annotated CDF files in the bioconductor	
		if(is.null(cdfPath))stop("This annotation needs a cdfPath argument") # but this can't really happen since we already checked for all being null
		if(!file.exists(cdfPath))stop(paste("The cdf file",cdfPath,"was not found"))
		if(file.info(cdfPath)[["isdir"]])stop(paste("The cdfPath",cdfPath,"is a directory"))
		if(is.null(annotation))stop("An annotation string is needed to locate the cdf and probe level files for CDF-style arrays") 
		annotation_type <- "hascdf"	
	}
	if(!annotation_type %in% c("hascdf","doesnothavecdf")){
		stop("Program could not determine if exon array or cdf-style arrays are investigated. Try giving a cdf_file or a pgf_file argument")
	}
	
	
	
	
	#if aptCelExtractPath is not given as argument
	if(is.null(aptCelExtractPath)){
		#..then check if aptCelExtractPath is in path
		aptCelExtractPath <- checkForFileInPath(c("apt-cel-extract","apt-cel-extract.exe"))
		if(length(aptCelExtractPath) > 0){
			aptCelExtractPath <- aptCelExtractPath[1] #only need one, in case of more hits
		}else{
			#...if not in path, check if it has been used before on this computer 
			locationFilePath <- file.path(.path.package("GeneRegionScan"),"configfiles","locations.txt")
			locationFile <- file(locationFilePath,"r")
			locationFileLines <- readLines(locationFile)
			close(locationFile)
			splitLocationFileLines <- strsplit(locationFileLines,"\t")
			for(line in splitLocationFileLines){
				if(line[1] == "aptCelExtractPath"){
					aptCelExtractPath <- line[3]
				}
			}
			
			#...if not used before, check if we can use the aptfiles in the exec folder of the package
			executableFileOverviewPath <- file.path(.path.package("GeneRegionScan"),"configfiles","aptextractfiles.txt")
			executableFileOverview <- read.table(executableFileOverviewPath,header=FALSE,sep="\t")
			machineType <- paste(Sys.info()[["sysname"]],Sys.info()[["machine"]])
			if(machineType %in% executableFileOverview[,1]){
				aptCelExtractName <- as.character(executableFileOverview[executableFileOverview[,1] == machineType,2])
				
				#Sometimes the file doesn't want to run on windows. We have to check that.
				if(Sys.info()[["sysname"]] == "Windows"){
					tempAptCelExtractPath <- file.path(.path.package("GeneRegionScan"),"exec",aptCelExtractName)
					ERRORLEVEL <- system(shQuote(tempAptCelExtractPath),intern=TRUE)
					if(length(ERRORLEVEL) > 20){ #if more than 20 lines -> safe to assume that the file worked and the help was printed
						aptCelExtractPath <- tempAptCelExtractPath
					}
				}else{
					aptCelExtractPath <- file.path(.path.package("GeneRegionScan"),"exec",aptCelExtractName)
				}
				
				
			}
			
			
			
			#...if still not found, we'll have to stop the show
			if(is.null(aptCelExtractPath)){
				stop("aptCelExtractPath was given as NULL, could not be found in path, and has not previously been located on this computer. Either make sure that the file apt-cel-extract or apt-cel-extract.exe is in path, or give the explicit path to it as an argument. The apt-cel-extract path can be downloaded from Affymetrix webpage")
			}
			
			#...if found, double check that it actually exists, is not a directory and is executable
			if(!file.exists(aptCelExtractPath)){
				stop(paste("The aptCelExtractPath was given as",aptCelExtractPath,"but no file was found here"))
			}
			if(file.info(aptCelExtractPath)[["isdir"]]){
				stop(paste("The aptCelExtractPath was given as",aptCelExtractPath,"but this is a directory"))
			}
			if(file.access(aptCelExtractPath,mode=1) != 0){
				stop(paste("The aptCelExtractPath was given as",aptCelExtractPath,"but this file is not executable"))
			}
		}
	}
	
	
	
	#saving the aptCelExtract, cdf, pgf and clf locations for next time
	#this functionality has been removed because it is not allowed to save settings between sessions in Bioconductor. 
	#it was re-inserted because other it seems that other packages do it as well
	if(file.access(locationFilePath,2) == 0){
		#locationFilePath <- file.path(.path.package("GeneRegionScan"),"configfiles","locations.txt")
		locationFile <- file(locationFilePath,"r")
		locationFileLines <- readLines(locationFile)
		close(locationFile)
		splitLocationFileLines <- strsplit(locationFileLines,"\t")
		newLocationFileLines <- vector()
		#delete previous lines from this annotation
		for(line in splitLocationFileLines){
			keepline <- TRUE
			if("aptCelExtractPath" == line[1])keepline <- FALSE
			if(annotation == line[2] & line[1] %in% c("cdfPath","clfPath","pgfPath"))keepline <- FALSE
			
			if(keepline){
				newLocationFileLines <- c(newLocationFileLines,paste(line[1],line[2],line[3],sep="\t"))
			}
		}	
		#add the new ones
		newLocationFileLines <- c(newLocationFileLines,paste("aptCelExtractPath","all",aptCelExtractPath,sep="\t"))
		if(!is.null(cdfPath))newLocationFileLines <- c(newLocationFileLines,paste("cdfPath",annotation,cdfPath,sep="\t"))
		if(!is.null(clfPath))newLocationFileLines <- c(newLocationFileLines,paste("clfPath",annotation,clfPath,sep="\t"))
		if(!is.null(pgfPath))newLocationFileLines <- c(newLocationFileLines,paste("pgfPath",annotation,pgfPath,sep="\t"))
		
		locationFile <- file(locationFilePath,"w")
		writeLines(newLocationFileLines,con=locationFile)
		close(locationFile)
		
	}
	if(file.access(locationFilePath,2) != 0 & verbose){
		print(paste("Would have saved the locations of annotation files for",annotation,"for next time, but did not have write permission for the configfile",locationFilePath))
	}
	
	
	#The following block is because I want users to both be able to submit vectors of probesets, and paths for files with probesets
	if(class(listOfProbesets[1]) == "numeric"|class(listOfProbesets[1]) == "integer")listOfProbesets <- as.character(listOfProbesets)
	if(length(listOfProbesets) == 1 & file.exists(listOfProbesets[1])){
		delete_probes_file <- FALSE
		if(!(readLines(listOfProbesets,n=2,warn=FALSE)[1] %in% c("probeset_id","probeset_name"))){
			oldlines <- readLines(listOfProbesets)
			probe_list_file <- file(listOfProbesets,"w")
			writeLines(c("probeset_id",oldlines),probe_list_file)
			close(probe_list_file)
		}
	}else{
		if(length(listOfProbesets) == 1){
			print("WARNING: the list of probesets only had one entry. This entry will be treated as a name of a probeset, but it is also possible that it is a missing file")
		}
		listOfProbesets <- c("probeset_id",listOfProbesets)
		probe_list_file <- file("probelist","w")
		writeLines(listOfProbesets,probe_list_file)
		close(probe_list_file)
		listOfProbesets <- paste(getwd(),.Platform[["file.sep"]],"probelist",sep="")
		delete_probes_file <- TRUE
	}
	
	vectorOfProbesets <- readLines(listOfProbesets)
	vectorOfProbesets <- vectorOfProbesets[2:length(vectorOfProbesets)]
	
	
	
	
	
	if(annotation_type == "hascdf"){
		annotation_line <- paste(" -d ",shQuote(cdfPath)," ",sep="")
		probelevel_package_name <- paste(annotation,"probe",sep="")
		cdf_package_name <- paste(annotation,"cdf",sep="")
		if(!require(probelevel_package_name,character.only=TRUE))stop(paste("The package",probelevel_package_name,"is required to obtain sequences of probes")) 
		if(!require(cdf_package_name,character.only=TRUE))stop(paste("The package",cdf_package_name,"is required to obtain dimensions of arrays"))
		dimensions <- as.list(get(paste(annotation,"dim",sep="")))
		if(verbose)print(paste("Now parsing the probesets received using",probelevel_package_name))
		probes <- subset(as.data.frame(get(probelevel_package_name)),as.data.frame(get(probelevel_package_name))[,"Probe.Set.Name"] %in% vectorOfProbesets)
		require(affy)
		get_probe_indices <- function(x){
			result <- xy2indices(as.numeric(x["x"]),as.numeric(x["y"]),cdf=cdf_package_name)
			return(result)
		}
		
		probenames <- apply(probes,1,get_probe_indices)
		
		probe_level_annotation <- cbind(probes[,"Probe.Set.Name" ],probes[,"sequence"])
		colnames(probe_level_annotation) <- c("probeset_name","sequence")
		rownames(probe_level_annotation) <- probenames
	}
	
	
	
	if(annotation_type == "doesnothavecdf"){
		
		
		annotation_line <- paste(" -c ",shQuote(clfPath)," -p ",shQuote(pgfPath)," ",sep="")
		if(verbose)print("Parsing the pgf file using readPgf. This might take a while.")
		probe_level_annotation <- getProbeLevelAnnotationForExonArrays(vectorOfProbesets,pgfPath=pgfPath)
	}
	
	if(sum(duplicated(rownames(probe_level_annotation))) > 0){
		duplicate_number<-sum(duplicated(rownames(probe_level_annotation)))
		stop(paste("There were",duplicate_number,"duplicate probe ids. The probe ids are not necessarily unique across the entire set, but should be unique in each probe set. Consider giving fewer probesets as argument."))
	}
	
	cel_files <- list.files(path=celfilePath,pattern="[.](c|C)(e|E)(l|L)$")
	
	cel_files_and_path <- paste(celfilePath,.Platform[["file.sep"]],cel_files,sep="")
	
	if(length(cel_files_and_path)<1){
		stop("Didn't find any cel files in specified celfilePath")
	}
	if(length(cel_files_and_path) == 1){
		stop("Only found one cel file in the specified celfilePath. This has not been implemented yet. Make a copy of the cel file if you really need to only read one")
	}
	
	
	if(verbose)print(paste("Found",length(cel_files_and_path),"cel files in celfilePath"))
	cel_list <- c("cel_files",cel_files_and_path)
	cel_list_file <- file("cellist","w")
	writeLines(cel_list,cel_list_file)
	close(cel_list_file)
	cel_list_path <- paste(getwd(),.Platform[["file.sep"]],"cellist",sep="")
	system(
			paste(
					shQuote(aptCelExtractPath),
					" -o ",
					shQuote(file.path(getwd(),"intensity_file.txt")),
					annotation_line,
					"  -a quant-norm,pm-only --probeset-ids ",
					shQuote(listOfProbesets),
					" --cel-files ",
					shQuote(cel_list_path),
					sep=""
			)
	)
	
	
	
	if(verbose)print("Now importing pgf-parsing data and intensity file")
	intensity_file <- paste(getwd(),.Platform[["file.sep"]],"intensity_file.txt",sep="")
	intensity <- try(read.table(intensity_file,header=TRUE,row.names=1,sep="\t"),silent=TRUE)
	
	if(class(intensity) == "try-error"){
		intensity <- try(read.table(intensity_file,header=TRUE,sep="\t"),silent=TRUE)
		if(class(intensity) == "try-error"){
			stop("There was an error reading the intensity file from apt. The cause is unknown")
		}else{
			duplicate_number<-sum(duplicated(intensity[,1]))
			stop(paste("Could not read the intensity file that was created by APT because of",duplicate_number,"duplicated row.names (some probes had same id). Consider entering fewer probesets to begin with."))
		}
	}
	
	
	intensity <- intensity[,7:ncol(intensity)]
	
	
	if(verbose){
		featurenames_before <- nrow(intensity)
		parsed_pgf_data_before <- nrow(probe_level_annotation)
	}
	
	intensity <- intensity[rownames(intensity) %in% rownames(probe_level_annotation),]
	probe_level_annotation <- probe_level_annotation[rownames(probe_level_annotation) %in% rownames(intensity),]
	
	
	if(verbose){
		featurenames_after <- nrow(intensity)
		parsed_pgf_data_after <- nrow(probe_level_annotation)
		if(featurenames_after/featurenames_before<1){
			percent_removed <- (1-round(featurenames_after/featurenames_before,3))*100
			if(percent_removed > 51 | percent_removed < 49){ #because we don't need to alarm people about removing MM-probes 
				print(paste("Removed",percent_removed,"% of the features returned by APT because they were not found when parsing the pgf file with readPgf"))
			}
		}
		if(parsed_pgf_data_after/parsed_pgf_data_before<1){
			percent_removed <- (1-round(parsed_pgf_data_after/parsed_pgf_data_before,3))*100
			print(paste("Removed",percent_removed,"% of the features returned by parsing the pgf file with readPgf because they were not found by APT"))
		}
	}
	
	title <- "ProbeLevelSet"	
	abstract <- "This ProbeLevelSet was generated using the getLocalProbeIntensities function of the bioconductor package geneRegionScan"
	
	experimentData  <-  new(
			"MIAME", 
			name = "", 
			lab = "", 
			contact = "", 
			title = title, 
			abstract = abstract, 
			url = "")
	
	
	expressionset <-  new(
			"ProbeLevelSet", 
			featureData = createFeatureData(probe_level_annotation), 
			exprs = data.matrix(intensity[rownames(probe_level_annotation),]), 
			experimentData = experimentData)
	
	# I wrapped this in a try block because I have observed some weird problems with 
	# accesing the annotation for different kinds of expressionsets and derivatives.
	if(class(try(expressionset@annotation)) == "character"){ 
		expressionset@annotation <- annotation
	}
	
	
	unlink(intensity_file)
	if(delete_probes_file)unlink(listOfProbesets)
	
	
	return(expressionset)
}





getProbesetsFromRegionOfInterest <- function(annotation,chromosome,start,end,pythonPath=NULL,transcriptClustersFile=NULL,mpsToPsFile=NULL){
	chromosome <- as.character(chromosome)
	
	if(is.null(mpsToPsFile) + is.null(transcriptClustersFile) == 1){
		stop("If supplying either mpsToPsFile or transcriptClustersFile you have to give them both")
	}
	
	if(is.null(mpsToPsFile) & is.null(transcriptClustersFile)){
		db_package_name <- paste(annotation,".db",sep="")
		if(db_package_name %in% .packages(all.available=TRUE)){ 
			require(AnnotationDbi)
			require(db_package_name,character.only=TRUE)
			chromosome_package <- paste(annotation,"CHR",sep="")
			if(substr(chromosome,1,3) == "chr")chromosome <- substr(chromosome,4,nchar(chromosome))
			probesets_on_correct_chromosome <- revmap(get(chromosome_package))[[chromosome]]
			chrloc_package <- paste(annotation,"CHRLOC",sep="")
			
			#yes it is slow, this following, but otherwise the database acts up
			probesets_of_interest <- vector()
			for(probeset in probesets_on_correct_chromosome){
				for(position in get(chrloc_package)[[probeset]]){
					found <- any(abs(position) > start & abs(position) < end)
					if(!is.na(found)){
						if(found){
							probesets_of_interest <- c(probesets_of_interest,probeset)
						}
					}
				}
			}
			
			
		}else{
			#received and annotation for which there was no db_package_name. Will check locations.txt if this
			#annotation has previously been run with mpsToPsFile and transcriptClustersFile.
			locationFilePath <- file.path(.path.package("GeneRegionScan"),"configfiles","locations.txt")
			if(file.access(locationFilePath,4) == 0){
				locationFile <- file(locationFilePath,"r")
				locationFileLines <- readLines(locationFile)
				close(locationFile)
				splitLocationFileLines <- strsplit(locationFileLines,"\t")
				for(line in splitLocationFileLines){
					if(annotation == line[2]){
						if("mpsToPsFile" == line[1])mpsToPsFile <- line[3]
						if("transcriptClustersFile" == line[1])transcriptClustersFile <- line[3]
					}
				}
				if(!is.null(transcriptClustersFile)){
					print(paste("transcriptClustersFile given as null, but the annotation",annotation,"has used this file before:",transcriptClustersFile))
				}
				if(!is.null(mpsToPsFile) ){
					print(paste("mpsToPsFile given as null, but the annotation",annotation,"has used this file before:",mpsToPsFile))
				}
				
			}
			if(is.null(transcriptClustersFile) | is.null(mpsToPsFile)){
				stop("The annotation was not found as any installed .db package. You either need to install the annotation.db package or specify a transcriptClustersFile and a mpsToPsFile in the case of exon arrays without .db packages")
			}
			
		}
	}
	if((!is.null(mpsToPsFile) & !is.null(transcriptClustersFile))){
		print("mpsToPsFile and transcriptClustersFile was detected. Switching to python parsing.")
		
		if(!file.exists(transcriptClustersFile))stop(paste("transcriptClustersFile file not found at",transcriptClustersFile))
		if(!file.exists(mpsToPsFile))stop(paste("mpsToPsFile file not found at",mpsToPsFile))
		
		if(is.null(pythonPath)){
			pythonPath <- checkForFileInPath(c("python","python.exe"))
			if(length(pythonPath) > 0){
				pythonPath <- pythonPath[1]
			}else{
				stop("Could not find python in path and it was not given as an argument. Do either")
			}
		}
		
		if(substr(chromosome,1,3) != "chr")chromosome <- paste("chr",chromosome,sep="")
		python_script <- c("def main(argv):",
				"    import os",
				"    transcriptClustersFile = argv[1]",
				"    mpsToPsFile = argv[2]",
				"    chr_name = argv[3]",
				"    start = int(argv[4])",
				"    end = int(argv[5])",
				"    transcript_clusters_file = open(transcriptClustersFile,'r')",
				"    transcript_clusters = transcript_clusters_file.readlines()",
				"    transcript_clusters_file.close()",
				"    transcript_clusters_of_interest = []",
				"    for transcript_cluster in transcript_clusters:",
				"        if transcript_cluster[0] is not '#':",
				"            split_cluster = transcript_cluster.split('\\\",\\\"')",
				"            if len(split_cluster) > 3:",
				"                if split_cluster[2] == chr_name:",
				"                    if (start < int(split_cluster[4]) and int(split_cluster[4]) < end) or (start < int(split_cluster[5]) and int(split_cluster[5]) < end):",
				"                        transcript_clusters_of_interest.append(split_cluster[0].lstrip('\\\"'))",
				"            else:",
				"                raise Exception('The transcriptClustersFile is not valid')",
				"    del transcript_clusters",
				"    mps_to_ps_file = open(mpsToPsFile,'r')",
				"    mps_to_ps = mps_to_ps_file.readlines()",
				"    mps_to_ps_file.close()",
				"    probesets_of_interest = []",
				"    for mps_to_ps_entry in mps_to_ps:",
				"        if mps_to_ps_entry[0] is not '#':",
				"            split_mps_to_ps_entry = mps_to_ps_entry.split('\t')",
				"            if split_mps_to_ps_entry[0] in transcript_clusters_of_interest:",
				"                probesets_of_interest_here = split_mps_to_ps_entry[2].split(' ')",
				"                probesets_of_interest = probesets_of_interest + probesets_of_interest_here",
				"    del mps_to_ps",
				"    probesets_of_interest = dict(map(lambda i: (i,1),probesets_of_interest)).keys()",
				"    output = open(os.path.join(os.getcwd(),'probesets_of_interest.txt'),'w')",
				"    output.write('probeset_id\\n')",
				"    for probeset_of_interest in probesets_of_interest:",
				"        if len(probeset_of_interest)>0:",
				"            output.write(probeset_of_interest+'\\n')",
				"    output.close()",
				"if __name__ == '__main__':",
				"    from sys import argv",
				"    main(argv)","","")
		
		python_file <- file(description="tempscript.py","w")
		writeLines(python_script,con=python_file)
		close(python_file)
		python_file_address <- file.path(getwd(),"tempscript.py")
		end <- as.integer(end)
		start <- as.integer(start)
		system(paste(shQuote(pythonPath),shQuote(python_file_address),shQuote(transcriptClustersFile),shQuote(mpsToPsFile),chromosome,start,end))
		if(!file.exists("probesets_of_interest.txt")){
			unlink(python_file_address)
			stop("R could not read the parsed data from Python. The cause of this error can probably be found in the error output a few lines above in the line beginning with Exception:")
		}
		probeset_file <- file("probesets_of_interest.txt","r")
		probesets_of_interest <- readLines(probeset_file)
		close(probeset_file)
		probesets_of_interest <- probesets_of_interest[2:length(probesets_of_interest)]
		unlink("probesets_of_interest.txt")
		unlink(python_file_address)
		
		
		
		#since it seems the read was succesfull with this annotation name and these mpsToPsFile and transcriptClustersFile, 
		#we will save the locations for future use
		#this functionality has been removed because it is not allowed to save settings between sessions in Bioconductor.
		#it was re-inserted because other it seems that other packages do it as well
		
		#saving the mpsToPsFile and transcriptClustersFile locations for next time 
		locationFilePath <- file.path(.path.package("GeneRegionScan"),"configfiles","locations.txt")
		if(file.access(locationFilePath,2) == 0){
			locationFile <- file(locationFilePath,"r")
			locationFileLines <- readLines(locationFile)
			close(locationFile)
			splitLocationFileLines <- strsplit(locationFileLines,"\t")
			newLocationFileLines <- vector()
			for(line in splitLocationFileLines){
				keepline <- TRUE
				if( !(annotation == line[2] & line[1] %in% c("mpsToPsFile","transcriptClustersFile")) ){
					newLocationFileLines <- c(newLocationFileLines,paste(line[1],line[2],line[3],sep="\t"))
				}
			}
			if(!is.null(mpsToPsFile))newLocationFileLines <- c(newLocationFileLines,paste("mpsToPsFile",annotation,mpsToPsFile,sep="\t"))
			if(!is.null(transcriptClustersFile))newLocationFileLines <- c(newLocationFileLines,paste("transcriptClustersFile",annotation,transcriptClustersFile,sep="\t"))
			
			locationFile <- file(locationFilePath,"w")
			writeLines(newLocationFileLines,con=locationFile)
			close(locationFile)
			
		}
		if(file.access(locationFilePath,2) != 0){
			print(paste("Would have saved the locations of mpsToPsFile and transcriptClustersFile files for",annotation,"for next time, but did not have write permission for the configfile",locationFilePath))
		}
		
		
		
	}
	
	return(probesets_of_interest)
}



getMetaprobesetsFromRegionOfInterest <- function(annotation,chromosome,start,end,pythonPath=NULL,transcriptClustersFile=NULL){
	chromosome <- as.character(chromosome)
	
	if(is.null(transcriptClustersFile)){
		
		#received an annotation for which there was no db_package_name. Will check locations.txt if this
		#annotation has previously been run with transcriptClustersFile.
		locationFilePath <- file.path(.path.package("GeneRegionScan"),"configfiles","locations.txt")
		if(file.access(locationFilePath,4) == 0){
			locationFile <- file(locationFilePath,"r")
			locationFileLines <- readLines(locationFile)
			close(locationFile)
			splitLocationFileLines <- strsplit(locationFileLines,"\t")
			for(line in splitLocationFileLines){
				if(annotation == line[2]){
					if("transcriptClustersFile" == line[1])transcriptClustersFile <- line[3]
				}
			}
			if(!is.null(transcriptClustersFile)){
				print(paste("transcriptClustersFile given as null, but the annotation",annotation,"has used this file before:",transcriptClustersFile))
			}
			
		}
		
		
	}
	
	if(is.null(transcriptClustersFile)){
		stop(paste("No transcriptClustersFile was given and nothing has been saved for",annotation,"in memory"))			
	}
	
	
	
	print("transcriptClustersFile was detected. Switching to python parsing.")
	
	if(!file.exists(transcriptClustersFile)){
		stop(paste("transcriptClustersFile file not found at",transcriptClustersFile))
	}
	
	
	if(is.null(pythonPath)){
		pythonPath <- checkForFileInPath(c("python","python.exe"))
		if(length(pythonPath) > 0){
			pythonPath <- pythonPath[1]
		}else{
			stop("Could not find python in path and it was not given as an argument. Do either")
		}
	}
	
	if(substr(chromosome,1,3) != "chr")chromosome <- paste("chr",chromosome,sep="")
	python_script <- c("def main(argv):",
			"    import os",
			"    transcriptClustersFile = argv[1]",
			"    chr_name = argv[2]",
			"    start = int(argv[3])",
			"    end = int(argv[4])",
			"    transcript_clusters_file = open(transcriptClustersFile,'r')",
			"    transcript_clusters = transcript_clusters_file.readlines()",
			"    transcript_clusters_file.close()",
			"    transcript_clusters_of_interest = []",
			"    for transcript_cluster in transcript_clusters:",
			"        if transcript_cluster[0] is not '#':",
			"            split_cluster = transcript_cluster.split('\\\",\\\"')",
			"            if len(split_cluster) > 3:",
			"                if split_cluster[2] == chr_name:",
			"                    if (start < int(split_cluster[4]) and int(split_cluster[4]) < end) or (start < int(split_cluster[5]) and int(split_cluster[5]) < end):",
			"                        transcript_clusters_of_interest.append(split_cluster[0].lstrip('\\\"'))",
			"            else:",
			"                raise Exception('The transcriptClustersFile is not valid')",
			"    del transcript_clusters",
			"    output = open(os.path.join(os.getcwd(),'metaprobesets_of_interest.txt'),'w')",
			"    output.write('probeset_id\\n')",
			"    for transcript_cluster_of_interest in transcript_clusters_of_interest:",
			"        if len(transcript_cluster_of_interest)>0:",
			"            output.write(transcript_cluster_of_interest+'\\n')",
			"    output.close()",
			"if __name__ == '__main__':",
			"    from sys import argv",
			"    main(argv)","","")
	
	python_file <- file(description="tempscript.py","w")
	writeLines(python_script,con=python_file)
	close(python_file)
	python_file_address <- file.path(getwd(),"tempscript.py")
	end <- as.integer(end)
	start <- as.integer(start)
	system(paste(shQuote(pythonPath),shQuote(python_file_address),shQuote(transcriptClustersFile),chromosome,start,end))
	if(!file.exists("metaprobesets_of_interest.txt")){
		unlink(python_file_address)
		stop("R could not read the parsed data from Python. The cause of this error can probably be found in the error output a few lines above in the line beginning with Exception:")
	}
	metaprobeset_file <- file("metaprobesets_of_interest.txt","r")
	metaprobesets_of_interest <- readLines(metaprobeset_file)
	close(metaprobeset_file)
	metaprobesets_of_interest <- metaprobesets_of_interest[2:length(metaprobesets_of_interest)]
	unlink("metaprobesets_of_interest.txt")
	unlink(python_file_address)
	
	
	
	#since it seems the read was succesfull with this annotation name and this transcriptClustersFile, 
	#we will save the locations for future use
	
	#saving the transcriptClustersFile locations for next time 
	#this functionality has been removed because it is not allowed to save settings between sessions in Bioconductor.
	#it was re-inserted because other it seems that other packages do it as well
	locationFilePath <- file.path(.path.package("GeneRegionScan"),"configfiles","locations.txt")
	if(file.access(locationFilePath,2) == 0){
		locationFile <- file(locationFilePath,"r")
		locationFileLines <- readLines(locationFile)
		close(locationFile)
		splitLocationFileLines <- strsplit(locationFileLines,"\t")
		newLocationFileLines <- vector()
		for(line in splitLocationFileLines){
			
			if( !(annotation == line[2] & line[1] %in% c("transcriptClustersFile")) ){
				newLocationFileLines <- c(newLocationFileLines,paste(line[1],line[2],line[3],sep="\t"))
			}
		}
		
		newLocationFileLines <- c(newLocationFileLines,paste("transcriptClustersFile",annotation,transcriptClustersFile,sep="\t"))
		
		locationFile <- file(locationFilePath,"w")
		writeLines(newLocationFileLines,con=locationFile)
		close(locationFile)
		
	}
	if(file.access(locationFilePath,2) != 0){
		print(paste("Would have saved the locations of transcriptClustersFile files for",annotation,"for next time, but did not have write permission for the configfile",locationFilePath))
	}
	
	
	
	return(metaprobesets_of_interest)
}





getProbeLevelAnnotationForExonArrays <- function(vectorOfProbesets,pgfPath){
	# When given a vectorOfProbesets and the location of a pgf file, this function will create
	# a dataframe with columns "probeset_name", "probe_name" and "sequence" with data for all probes
	# found in the vectorOfProbesets. If any probesets are given that are not found in the pgf file
	# the function will stop.
	#
	
	
	
	parsed_pgf <- readPgf(pgfPath) #this can produce a 
	
	probeset_indices_logical <- parsed_pgf[["probesetId"]] %in% vectorOfProbesets
	
	probeset_indices <- which(probeset_indices_logical)
	
	
	#Might as well stop the function if it doesn't find all probesets, because APT will crash later in any case
	if(sum(probeset_indices_logical)<length(vectorOfProbesets)){
		names_of_missing_probesets <- vectorOfProbesets[!vectorOfProbesets %in% parsed_pgf[["probesetId"]][probeset_indices_logical]]
		number_of_missing_probesets <- length(names_of_missing_probesets)
		if(number_of_missing_probesets > 50){
			description_of_missing_probesets <- paste("The first 50 are:",paste(names_of_missing_probesets[1:50],collapse=", "))
		}else{
			description_of_missing_probesets <- paste("They are:",paste(names_of_missing_probesets,collapse=", "))
		}
		stop(paste(number_of_missing_probesets,"probesets from the vectorOfProbesets, was not found in the pgf file.",description_of_missing_probesets))
	}
	
	probe_indices <- vector()
	probeset_name <- vector()
	
	for(probeset_index in probeset_indices){
		atom_indices <- parsed_pgf[["probesetStartAtom"]][probeset_index]:(parsed_pgf[["probesetStartAtom"]][probeset_index+1]-1)
		for(atom_index in atom_indices){
			probe_indices_here <- parsed_pgf[["atomStartProbe"]][atom_index]:(parsed_pgf[["atomStartProbe"]][atom_index+1]-1)
			probe_indices <- c(probe_indices,probe_indices_here)
			probeset_name <- c(probeset_name,rep(parsed_pgf[["probesetId"]][probeset_index],length(atom_index)))
		}
	}
	
	
	probe_name <- parsed_pgf[["probeId"]][probe_indices]
	sequence <- parsed_pgf[["probeSequence"]][probe_indices]
	
	result <- cbind(probeset_name,probe_name,sequence)
	colnames(result) <- c("probeset_name","probe_name","sequence")
	rownames(result) <- result[,"probe_name"]
	return(result)
}











checkForFileInPath <- function(filenames){
	#small function that takes a vector of filenames, checks if any is in path, and returns the full path of thoose that are.
	# if nothing is found, it will return NULL
	return_value <- vector()
	for(filename in filenames){
		pathdirs <- strsplit(Sys.getenv()[["PATH"]],.Platform[["path.sep"]])[[1]]
		for(pathdir in pathdirs){
			if(!is.na(file.info(pathdir)$isdir)){ #this sometimes happens on some linux computers. Weird, but the warning should be avoided.
				if(filename %in% list.files(pathdir)){
					return_value <- c(return_value,paste(pathdir,filename,sep=.Platform$file.sep))
				}
			}
		}
	}
	if(length(return_value) == 0)return_value <- NULL
	return(return_value)
}









getLocalMetaprobeIntensities <- function(celfilePath,analysis="rma",metaProbeSetsFile=NULL,annotation=NULL,aptProbesetSummarizePath=NULL,pgfPath=NULL,clfPath=NULL,cdfPath=NULL,verbose=TRUE){
	if(is.null(cdfPath) & is.null(pgfPath) & is.null(clfPath) & is.null(annotation)){
		stop("You need to provide at least one of the following: pgfPath, clfPath, cdfPath or annotation - and probably more, but error messages will tell you which")
	}
	if(!is.null(cdfPath) & (!is.null(pgfPath) | !is.null(clfPath))){
		stop("Do not specify both cdfPath and a clfPath or pgfPath. pgfPath and clfPath is for exon arrays that do not use cdf files")
	}
	
	if(!file.exists(celfilePath))stop(paste("The celfilePath",celfilePath,"was not found"))
	if(!file.info(celfilePath)[["isdir"]])stop(paste("The celfilePath",celfilePath,"was not identified as a directory. You must specify a directory."))
	if(substr(celfilePath,nchar(celfilePath),nchar(celfilePath)) == .Platform[["file.sep"]]){
		celfilePath <- substr(celfilePath,1,nchar(celfilePath)-length(.Platform[["file.sep"]]))
	}
	celfilePath <- path.expand(celfilePath)
	
	if(!is.null(annotation)){
		if(length(grep("pgfPath",annotation)) + length(grep("clfPath",annotation)) + length(grep("cdfPath",annotation)) > 0){
			stop("For some obscure reason the words cdfPath, pgfPath or clfPath appears in the annotation name. This will cause the program to crash later. Please change the annotation name.")
		}
	}
	
	#if all cdf, pgf and clf is set to null, but annotation is not, we'll assume it is a lazy user
	#and see if this type of annotation has been used before.
	locationFilePath <- file.path(.path.package("GeneRegionScan"),"configfiles","locations.txt")
	if(is.null(cdfPath) & is.null(pgfPath) & is.null(clfPath) & !is.null(annotation) & (file.access(locationFilePath,4) == 0)){
		locationFile <- file(locationFilePath,"r")
		locationFileLines <- readLines(locationFile)
		close(locationFile)
		splitLocationFileLines <- strsplit(locationFileLines,"\t")
		for(line in splitLocationFileLines){
			if(annotation == line[2]){
				if("pgfPath" == line[1])pgfPath <- line[3]
				if("clfPath" == line[1])clfPath <- line[3]
				if("cdfPath" == line[1])cdfPath <- line[3]
			}
		}
		if(!is.null(cdfPath) & verbose)print(paste("cdfPath given as null, but the annotation",annotation,"has used this file before:",cdfPath))
		if(!is.null(clfPath) & verbose)print(paste("clfPath given as null, but the annotation",annotation,"has used this file before:",clfPath))
		if(!is.null(pgfPath) & verbose)print(paste("pgfPath given as null, but the annotation",annotation,"has used this file before:",pgfPath))
	}
	#
	annotation_type <- ""
	if(is.null(cdfPath) & is.null(pgfPath) & is.null(clfPath))stop(paste("For the annotation",annotation,"you need to provide at least one of the following: pgfPath, clfPath, cdfPath"))
	if(is.null(cdfPath)){
		#assuming exon type array
		if(is.null(pgfPath)|is.null(clfPath))stop("This annotation needs both a pgfPath argument and a cdfPath argument")
		if(!file.exists(pgfPath)|!file.exists(clfPath))stop(paste("The pgf or clf file was not found"))
		if(file.info(pgfPath)[["isdir"]]|file.info(clfPath)[["isdir"]])stop(paste("The pgf or clf specified is a directory, not a file"))
		annotation_type <- "doesnothavecdf"
	}
	if(is.null(clfPath)|is.null(pgfPath)){
		#assuming 3'IVT type array (ie with already annotated CDF files in the bioconductor	
		if(is.null(cdfPath))stop("This annotation needs a cdfPath argument") # but this can't really happen since we already checked for all being null
		if(!file.exists(cdfPath))stop(paste("The cdf file",cdfPath,"was not found"))
		if(file.info(cdfPath)[["isdir"]])stop(paste("The cdfPath",cdfPath,"is a directory"))
		if(is.null(annotation))stop("An annotation string is needed to locate the cdf and probe level files for CDF-style arrays") 
		annotation_type <- "hascdf"	
	}
	if(!annotation_type %in% c("hascdf","doesnothavecdf"))stop("Program could not determine if exon array or cdf-style arrays are investigated. Try giving a cdf_file or a pgf_file argument")
	
	
	
	#checking the metaProbeSetsFile
	if(!is.null(metaProbeSetsFile)){
		if(!file.exists(metaProbeSetsFile)){
			stop(paste("The metaProbeSetsFile was given as",metaProbeSetsFile,"but no file was found here"))
		}
		if(file.info(metaProbeSetsFile)[["isdir"]]){
			stop(paste("The metaProbeSetsFile was given as",metaProbeSetsFile,"but this is a directory"))
		}
		
		metaProbeSetline<-paste(" -m ",shQuote(metaProbeSetsFile)," ",sep="")
	}else{
		metaProbeSetline<-""
	}
	
	
	#if aptProbesetSummarizePath is not given as argument
	if(is.null(aptProbesetSummarizePath)){
		#..then check if aptProbesetSummarizePath is in path
		aptProbesetSummarizePath <- checkForFileInPath(c("apt-probeset-summarize","apt-probeset-summarize.exe"))
		if(length(aptProbesetSummarizePath) > 0){
			aptProbesetSummarizePath <- aptProbesetSummarizePath[1] #only need one, in case of more hits
		}else{
			#...if not in path, check if it has been used before on this computer 
			locationFilePath <- file.path(.path.package("GeneRegionScan"),"configfiles","locations.txt")
			locationFile <- file(locationFilePath,"r")
			locationFileLines <- readLines(locationFile)
			close(locationFile)
			splitLocationFileLines <- strsplit(locationFileLines,"\t")
			for(line in splitLocationFileLines){
				if(line[1] == "aptProbesetSummarizePath"){
					aptProbesetSummarizePath <- line[3]
				}
			}
			
			
			#...if not used before, check if we can use the aptfiles in the exec folder of the package
			executableFileOverviewPath <- file.path(.path.package("GeneRegionScan"),"configfiles","aptsummarizefiles.txt")
			executableFileOverview <- read.table(executableFileOverviewPath,header=FALSE,sep="\t")
			machineType <- paste(Sys.info()[["sysname"]],Sys.info()[["machine"]])
			if(machineType %in% executableFileOverview[,1]){
				aptCelSummarizetName <- as.character(executableFileOverview[executableFileOverview[,1] == machineType,2])
				
				#Sometimes the file doesn't want to run on windows. We have to check that.
				if(Sys.info()[["sysname"]] == "Windows"){
					tempAptCelSummarizePath <- file.path(.path.package("GeneRegionScan"),"exec",aptCelSummarizetName)
					ERRORLEVEL <- system(shQuote(tempAptCelSummarizePath),intern=TRUE)
					if(length(ERRORLEVEL) > 20){ #if more than 20 lines -> safe to assume that the file worked and the help was printed
						aptProbesetSummarizePath <- tempAptCelSummarizePath
					}
				}else{
					aptProbesetSummarizePath <- file.path(.path.package("GeneRegionScan"),"exec",aptCelSummarizetName)
				}
				
				
			}
			
			#...if still not found, we'll have to stop the show
			if(is.null(aptProbesetSummarizePath)){
				stop("aptProbesetSummarizePath was given as NULL, could not be found in path, and has not previously been located on this computer. Either make sure that the file apt-probeset-summarize or apt-probeset-summarize.exe is in path, or give the explicit path to it as an argument. The apt-probeset-summarize can be downloaded from Affymetrix webpage")
			}
			
			#...if found, double check that it actually exists, is not a directory and is executable
			if(!file.exists(aptProbesetSummarizePath)){
				stop(paste("The aptProbesetSummarizePath was given as",aptProbesetSummarizePath,"but no file was found here"))
			}
			if(file.info(aptProbesetSummarizePath)[["isdir"]]){
				stop(paste("The aptProbesetSummarizePath was given as",aptProbesetSummarizePath,"but this is a directory"))
			}
			if(file.access(aptProbesetSummarizePath,mode=1) != 0){
				stop(paste("The aptProbesetSummarizePath was given as",aptProbesetSummarizePath,"but this file is not executable"))
			}
		}
	}
	
	
	
	#saving the aptCelExtract, cdf, pgf and clf locations for next time 
	#this functionality has been removed because it is not allowed to save settings between sessions in Bioconductor.
	#it was re-inserted because other it seems that other packages do it as well
	if(file.access(locationFilePath,2) == 0){
		#locationFilePath <- file.path(.path.package("GeneRegionScan"),"configfiles","locations.txt")
		locationFile <- file(locationFilePath,"r")
		locationFileLines <- readLines(locationFile)
		close(locationFile)
		splitLocationFileLines <- strsplit(locationFileLines,"\t")
		newLocationFileLines <- vector()
		#delete previous lines from this annotation
		for(line in splitLocationFileLines){
			keepline <- TRUE
			if("aptProbesetSummarizePath" == line[1])keepline <- FALSE
			if(annotation == line[2] & line[1] %in% c("cdfPath","clfPath","pgfPath"))keepline <- FALSE
			
			if(keepline){
				newLocationFileLines <- c(newLocationFileLines,paste(line[1],line[2],line[3],sep="\t"))
			}
		}	
		#add the new ones
		newLocationFileLines <- c(newLocationFileLines,paste("aptProbesetSummarizePath","all",aptProbesetSummarizePath,sep="\t"))
		if(!is.null(cdfPath))newLocationFileLines <- c(newLocationFileLines,paste("cdfPath",annotation,cdfPath,sep="\t"))
		if(!is.null(clfPath))newLocationFileLines <- c(newLocationFileLines,paste("clfPath",annotation,clfPath,sep="\t"))
		if(!is.null(pgfPath))newLocationFileLines <- c(newLocationFileLines,paste("pgfPath",annotation,pgfPath,sep="\t"))
		
		locationFile <- file(locationFilePath,"w")
		writeLines(newLocationFileLines,con=locationFile)
		close(locationFile)
		
	}
	if(file.access(locationFilePath,2) != 0 & verbose){
		print(paste("Would have saved the locations of annotation files for",annotation,"for next time, but did not have write permission for the configfile",locationFilePath))
	}
	
	
	
	
	if(annotation_type == "hascdf"){
		annotation_line <- paste(" -d ",shQuote(cdfPath)," ",sep="")
	}
	
	if(annotation_type == "doesnothavecdf"){
		annotation_line <- paste(" -c ",shQuote(clfPath)," -p ",shQuote(pgfPath)," ",sep="")
	}
	
	
	
	
	cel_files <- list.files(path=celfilePath,pattern="[.](c|C)(e|E)(l|L)$")
	
	cel_files_and_path <- paste(celfilePath,.Platform[["file.sep"]],cel_files,sep="")
	
	if(length(cel_files_and_path)<1){
		stop("Didn't find any cel files in specified celfilePath")
	}
	if(length(cel_files_and_path) == 1){
		stop("Only found one cel file in the specified celfilePath. This has not been implemented yet. Make a copy of the cel file if you really need to only read one")
	}
	
	
	if(verbose)print(paste("Found",length(cel_files_and_path),"cel files in celfilePath"))
	cel_list <- c("cel_files",cel_files_and_path)
	cel_list_file <- file("cellist","w")
	writeLines(cel_list,cel_list_file)
	close(cel_list_file)
	cel_list_path <- paste(getwd(),.Platform[["file.sep"]],"cellist",sep="")
	results_file_path<-file.path(getwd(),"temp_intensity_files")
	
	system(
			paste(
					shQuote(aptProbesetSummarizePath),
					" -a ",
					analysis,
					annotation_line,					
					" -o ",
					shQuote(results_file_path),
					metaProbeSetline,
					" --cel-files ",
					shQuote(cel_list_path),
					sep=""
			)
	)
	
	
	if(verbose)print("APT analysis finished - now preparing values into bioconductor expressionset") 
	
	results_file_name_and_path<-file.path(results_file_path,paste(analysis,".summary.txt",sep=""))
	results_file<-file(results_file_name_and_path,open="r")
	results<-readLines(results_file)
	close(results_file)
	
	
	rma_report<-vector()
	for(i in 1:length(results)){
		line<-results[i]
		if(substr(line,1,1)=="#"){
			rma_report<-c(rma_report,line)
		}else{
			break
		}
	}
	skip_these<-i-1
	exprs<-as.matrix(read.table(results_file_name_and_path,sep="\t",skip=skip_these,header=TRUE,row.names=1))
	
	summary_file<-file.path(results_file_path,"apt-probeset-summarize.log")
	if(file.exists(summary_file)){
		summarize_log_file<-file(summary_file,open="r")
		summarize_log<-readLines(summarize_log_file)
		close(summarize_log_file)
		notes<-list(summarize_log,rma_report)
		names(notes)<-c("summarize_log","rma_report")
	}else{
		notes<-list(rma_report)
		names(notes)<-c("rma_report")
		
	}
	
	title <- "ExpressionSet"
	abstract <- paste("This ExpressionSet was generated using the getLocalMetaprobeIntensities function of the bioconductor package geneRegionScan, using",analysis,"analysis directly from cel files in",celfilePath,"on",date())
	experimentData  <-  new(
			"MIAME", 
			name = "", 
			lab = "", 
			contact = "", 
			title = title, 
			abstract = abstract, 
			url = "",
			other = notes)
	
	expressionset <- new("ExpressionSet", exprs = exprs, experimentData = experimentData, annotation = annotation)
	
	if(verbose)print("Deleting the files: apt-probeset-summarize.log, cel_list.txt, rma.summary.txt, and rma.report.txt - since they are now in the expressionset")
	file.remove(paste(results_file_path,c("apt-probeset-summarize.log",paste(analysis,".summary.txt",sep=""),paste(analysis,".report.txt",sep="")),sep=.Platform$file.sep))
	
	file.remove(cel_list_path)
	return(expressionset)
}




getProbesetsFromMetaprobeset <- function(annotation,metaprobesets,pythonPath=NULL,mpsToPsFile=NULL){
	
	
	if(is.null(mpsToPsFile)){
		#Will check locations.txt if this
		#annotation has previously been run with mpsToPsFile and transcriptClustersFile.
		locationFilePath <- file.path(.path.package("GeneRegionScan"),"configfiles","locations.txt")
		if(file.access(locationFilePath,4) == 0){
			locationFile <- file(locationFilePath,"r")
			locationFileLines <- readLines(locationFile)
			close(locationFile)
			splitLocationFileLines <- strsplit(locationFileLines,"\t")
			for(line in splitLocationFileLines){
				if(annotation == line[2]){
					if("mpsToPsFile" == line[1])mpsToPsFile <- line[3]
				}
			}
			if(!is.null(mpsToPsFile) ){
				print(paste("mpsToPsFile given as null, but the annotation",annotation,"has used this file before:",mpsToPsFile))
			}
			
		}
		if(is.null(mpsToPsFile)){
			stop("You need to specify a mpsToPsFile")
		}
		
	}
	
	if(!is.null(mpsToPsFile) ){
		
		
		if(!file.exists(mpsToPsFile))stop(paste("mpsToPsFile file not found at",mpsToPsFile))
		
		if(is.null(pythonPath)){
			pythonPath <- checkForFileInPath(c("python","python.exe"))
			if(length(pythonPath) > 0){
				pythonPath <- pythonPath[1]
			}else{
				stop("Could not find python in path and it was not given as an argument. Do either")
			}
		}
		
		
		python_script <- c("def main(argv):",
				"    import os",
				"    mpsToPsFile = argv[1]",
				"    metaprobesets = argv[2:(len(argv))]",
				"    mps_to_ps_file = open(mpsToPsFile,'r')",
				"    mps_to_ps = mps_to_ps_file.readlines()",
				"    mps_to_ps_file.close()",
				"    probesets_of_interest = []",
				"    for mps_to_ps_entry in mps_to_ps:",
				"        if mps_to_ps_entry[0] is not '#':",
				"            split_mps_to_ps_entry = mps_to_ps_entry.split('\t')",
				"            if split_mps_to_ps_entry[0] in metaprobesets:",
				"                probesets_of_interest_here = split_mps_to_ps_entry[2].split(' ')",
				"                probesets_of_interest = probesets_of_interest + probesets_of_interest_here",
				"    del mps_to_ps",
				"    probesets_of_interest = dict(map(lambda i: (i,1),probesets_of_interest)).keys()",
				"    output = open(os.path.join(os.getcwd(),'probesets_of_interest.txt'),'w')",
				"    output.write('probeset_id\\n')",
				"    for probeset_of_interest in probesets_of_interest:",
				"        if len(probeset_of_interest)>0:",
				"            output.write(probeset_of_interest+'\\n')",
				"    output.close()",
				"if __name__ == '__main__':",
				"    from sys import argv",
				"    main(argv)","","")
		
		python_file <- file(description="tempscript.py","w")
		writeLines(python_script,con=python_file)
		close(python_file)
		python_file_address <- file.path(getwd(),"tempscript.py")
		system(paste(shQuote(pythonPath),shQuote(python_file_address),shQuote(mpsToPsFile),paste(metaprobesets,collapse=" ")))
		if(!file.exists("probesets_of_interest.txt")){
			unlink(python_file_address)
			stop("R could not read the parsed data from Python. The cause of this error can probably be found in the error output a few lines above in the line beginning with Exception:")
		}
		probeset_file <- file("probesets_of_interest.txt","r")
		probesets_of_interest <- readLines(probeset_file)
		close(probeset_file)
		probesets_of_interest <- probesets_of_interest[2:length(probesets_of_interest)]
		unlink("probesets_of_interest.txt")
		unlink(python_file_address)
		
		
		
		#since it seems the read was succesfull with this annotation name and these mpsToPsFile, 
		#we will save the locations for future use
		#this functionality has been removed because it is not allowed to save settings between sessions in Bioconductor.
		#it was re-inserted because other it seems that other packages do it as well
		
		#saving the mpsToPsFilelocations for next time 
		locationFilePath <- file.path(.path.package("GeneRegionScan"),"configfiles","locations.txt")
		if(file.access(locationFilePath,2) == 0){
			locationFile <- file(locationFilePath,"r")
			locationFileLines <- readLines(locationFile)
			close(locationFile)
			splitLocationFileLines <- strsplit(locationFileLines,"\t")
			newLocationFileLines <- vector()
			for(line in splitLocationFileLines){
				if( !(annotation == line[2] & line[1] == "mpsToPsFile") ){
					newLocationFileLines <- c(newLocationFileLines,paste(line[1],line[2],line[3],sep="\t"))
				}
			}
			if(!is.null(mpsToPsFile))newLocationFileLines <- c(newLocationFileLines,paste("mpsToPsFile",annotation,mpsToPsFile,sep="\t"))
			
			
			locationFile <- file(locationFilePath,"w")
			writeLines(newLocationFileLines,con=locationFile)
			close(locationFile)
			
		}
		if(file.access(locationFilePath,2) != 0){
			print(paste("Would have saved the locations of mpsToPsFile and transcriptClustersFile files for",annotation,"for next time, but did not have write permission for the configfile",locationFilePath))
		}
		
		
		
	}
	
	return(probesets_of_interest)
}






readFASTA_replacement<-function(path){
	#small replacement function for the deprecated readFASTA
	DNAStringSet<-read.DNAStringSet(path)
	readFASTAformatList<-list()
	for(i in 1:length(DNAStringSet)){
		readFASTAformatList[[i]]<-list()
		
		readFASTAformatList[[i]][["desc"]]<- names(as.character(DNAStringSet))[i]
		readFASTAformatList[[i]][["seq"]]<- as.character(as.character(DNAStringSet)[i])
	}
	return(readFASTAformatList)
}
