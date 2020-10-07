#' @name geneNames
#' @inherit AcidGenerics::geneNames
#' @note Updated 2019-07-28.
#'
#' @inheritParams AcidRoxygen::params
#' @param ... Additional arguments.
#'
#' @examples
#' data(RangedSummarizedExperiment, package = "AcidTest")
#' rse <- RangedSummarizedExperiment
#'
#' ## SummarizedExperiment ====
#' head(geneNames(rse))
NULL



## Updated 2019-07-22.
`geneNames,Vector` <-  # nolint
    function(object) {
        Gene2Symbol(object, format = "makeUnique")[["geneName"]]
    }



## Updated 2019-07-22.
`geneNames,GRanges` <-  # nolint
    `geneNames,Vector`



#' @rdname geneNames
#' @export
setMethod(
    f = "geneNames",
    signature = signature("GRanges"),
    definition = `geneNames,GRanges`
)



## Updated 2019-07-22.
`geneNames,SummarizedExperiment` <-  # nolint
    `geneNames,Vector`



#' @rdname geneNames
#' @export
setMethod(
    f = "geneNames",
    signature = signature("SummarizedExperiment"),
    definition = `geneNames,SummarizedExperiment`
)
