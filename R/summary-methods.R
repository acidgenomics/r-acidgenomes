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
#' ## Gene2Symbol ====
#' df <- DataFrame(
#'     "geneId" = c(
#'         "ENSG00000228572.7",
#'         "ENSG00000182378.14"
#'     ),
#'     "geneName" = c(
#'         "AL954722.1",
#'         "PLCXD1"
#'     )
#' )
#' g2s <- Gene2Symbol(df)
#' summary(g2s)
#'
#' ## Tx2Gene ====
#' df <- DataFrame(
#'     "txId" = c(
#'         "ENST00000635602.1",
#'         "ENST00000635506.1"
#'     ),
#'     "geneId" = c(
#'         "ENSG00000283061.1",
#'         "ENSG00000283061.1"
#'     )
#' )
#' t2g <- Tx2Gene(df)
#' summary(t2g)
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
