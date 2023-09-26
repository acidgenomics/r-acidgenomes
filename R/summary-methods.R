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
#' ## GeneToSymbol ====
#' df <- S4Vectors::DataFrame(
#'     "geneId" = c(
#'         "ENSG00000228572.7",
#'         "ENSG00000182378.14"
#'     ),
#'     "geneName" = c(
#'         "AL954722.1",
#'         "PLCXD1"
#'     )
#' )
#' g2s <- GeneToSymbol(df)
#' summary(g2s)
#'
#' ## TxToGene ====
#' df <- S4Vectors::DataFrame(
#'     "txId" = c(
#'         "ENST00000635602.1",
#'         "ENST00000635506.1"
#'     ),
#'     "geneId" = c(
#'         "ENSG00000283061.1",
#'         "ENSG00000283061.1"
#'     )
#' )
#' t2g <- TxToGene(df)
#' summary(t2g)
NULL



## Updated 2021-03-03.
`summary,GeneToSymbol` <- # nolint
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
            "date" = m[["date"]]
        ))
    }



## Updated 2021-03-03.
`summary,TxToGene` <- # nolint
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
            "date" = m[["date"]]
        ))
    }



#' @rdname summary
#' @export
setMethod(
    f = "summary",
    signature = signature(object = "GeneToSymbol"),
    definition = `summary,GeneToSymbol`
)

#' @rdname summary
#' @export
setMethod(
    f = "summary",
    signature = signature(object = "TxToGene"),
    definition = `summary,TxToGene`
)
