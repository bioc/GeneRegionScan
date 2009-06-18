setClass("ProbeLevelSet", contains="ExpressionSet")

setValidity("ProbeLevelSet", function(object) {
    assayDataValidMembers(assayData(object))
})

