#' Object summaries
#'
#' @name summary
#' @keywords internal
#' @note Updated 2019-11-19.
#'
#' @inheritParams AcidRoxygen::params
#' @param ... Additional arguments.
#'
#' @return Invisible `NULL`.
#'
#' @examples
#' data(
#'     RangedSummarizedExperiment,
#'     SummarizedExperiment_transcripts,
#'     package = "AcidTest"
#' )
#'
#' ## Gene2Symbol ====
#' object <- Gene2Symbol(RangedSummarizedExperiment)
#' summary(object)
#'
#' ## Tx2Gene ====
#' object <- Tx2Gene(SummarizedExperiment_transcripts)
#' summary(object)
NULL



## Updated 2021-02-01.
`summary,Gene2Symbol` <-  # nolint
    function(object) {
        m <- metadata(object)
        showSlotInfo(list(
            "genes" = length(unique(object[["geneId"]])),
            "symbols" = length(unique(object[["geneName"]])),
            "format" = m[["format"]],
            "organism" = m[["organism"]],
            "provider" = m[["provider"]],
            "genomeBuild" = m[["genomeBuild"]],
            "release" = m[["release"]],
            "annotationHub" = m[["annotationHubId"]],
            "acidGenomes" = as.character(m[["acidGenomes"]]),
            "date" = m[["date"]]
        ))
    }



#' @rdname summary
#' @export
setMethod(
    f = "summary",
    signature = signature("Gene2Symbol"),
    definition = `summary,Gene2Symbol`
)



## Updated 2021-02-01.
`summary,Tx2Gene` <-  # nolint
    function(object) {
        m <- metadata(object)
        showSlotInfo(list(
            "transcripts" = length(unique(object[["txId"]])),
            "genes" = length(unique(object[["geneId"]])),
            "organism" = m[["organism"]],
            "provider" = m[["provider"]],
            "genomeBuild" = m[["genomeBuild"]],
            "release" = m[["ensemblRelease"]],
            "annotationHub" = m[["id"]],
            "acidGenomes" = as.character(m[["acidGenomes"]]),
            "date" = m[["date"]]
        ))
    }



#' @rdname summary
#' @export
setMethod(
    f = "summary",
    signature = signature("Tx2Gene"),
    definition = `summary,Tx2Gene`
)
