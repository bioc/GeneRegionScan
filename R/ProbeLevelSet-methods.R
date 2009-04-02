# returns the sequence of the specified probes.
setGeneric("getSequence", function(object, id) standardGeneric("getSequence"))
setMethod("getSequence", "ProbeLevelSet",
function(object, id) {
	if(missing(id)) {
		result<-pData(featureData(object))[, "sequence"]
		names(result)<-rownames(pData(featureData(object)))
	} else {
		result<-pData(featureData(object))[id, "sequence"]
		names(result)<-id
	}
	return(result)
})

createFeatureData <- function(x) {
	featuredata <- data.frame(probeset_name = x[, "probeset_name"], sequence = x[, "sequence"],stringsAsFactors=FALSE)
	rownames(featuredata) <- rownames(x)
    featuredata.feMet <- data.frame(labelDescription = c("probeset id", "probe sequence"), row.names = c("probesetId", "probeSeq"))
    new("AnnotatedDataFrame", data = featuredata, varMetadata = featuredata.feMet)
}
