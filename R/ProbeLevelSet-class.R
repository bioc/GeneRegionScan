setClass("ProbeLevelSet", contains="ExpressionSet")

#
setMethod("initialize", "ProbeLevelSet",
function(.Object,
    #phenoData = new("AnnotatedDataFrame"),
    #featureData = new("AnnotatedDataFrame"),
    #experimentData = new("MIAME"),
    #annotation = character(),
    ...
) {
    callNextMethod(.Object,
        #phenoData = phenoData,
        #featureData = featureData,
        #experimentData = experimentData,
        #annotation = annotation,
        ...
    )
})

setValidity("ProbeLevelSet", function(object) {
    assayDataValidMembers(assayData(object))
})

