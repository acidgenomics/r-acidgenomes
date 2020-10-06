#' @name stripTranscriptVersions
#' @inherit AcidGenerics::stripTranscriptVersions
#' @note Updated 2019-10-09.
#'
#' @inheritParams AcidRoxygen::params
#' @param ... Additional arguments.
#'
#' @note This method is strict, and will only strip Ensembl transcript IDs
#'   beginning with "ENS".
#'
#' @examples
#' ## Ensembl (modify; contains versions)
#' stripTranscriptVersions(c(
#'     "ENSMUST00000000001.1",
#'     "ENSMUST00000000001-1",
#'     "ENSMUST00000000001_1"
#' ))
#'
#' ## WormBase (keep; doesn't contain versions)
#' stripTranscriptVersions("cTel79B.1")
NULL



#' @rdname stripTranscriptVersions
#' @name stripTranscriptVersions
#' @importFrom AcidGenerics stripTranscriptVersions
#' @usage stripTranscriptVersions(object, ...)
#' @export
NULL



## Pattern matching against Ensembl transcript (and gene) IDs.
##
## Example prefixes: ENST (human); ENSMUST (mouse).
## `:punct:` will match `-` or `_` here.
##
## See also:
## - http://www.ensembl.org/info/genome/stable_ids/index.html
##
## Updated 2019-10-07.
`stripTranscriptVersions,character` <-  # nolint
    function(object) {
        assert(isCharacter(object))
        gsub(
            pattern = "^(ENS.*[GT][[:digit:]]{11})[[:punct:]][[:digit:]]+$",
            replacement = "\\1",
            x = object
        )
    }



#' @rdname stripTranscriptVersions
#' @export
setMethod(
    f = "stripTranscriptVersions",
    signature = signature("character"),
    definition = `stripTranscriptVersions,character`
)



## Updated 2019-07-22.
`stripTranscriptVersions,matrix` <-  # nolint
    function(object) {
        rownames <- rownames(object)
        rownames <- stripTranscriptVersions(rownames)
        rownames(object) <- rownames
        object
    }



#' @rdname stripTranscriptVersions
#' @export
setMethod(
    f = "stripTranscriptVersions",
    signature = signature("matrix"),
    definition = `stripTranscriptVersions,matrix`
)



## Updated 2020-01-30.
`stripTranscriptVersions,Matrix` <-  # nolint
    `stripTranscriptVersions,matrix`



#' @rdname stripTranscriptVersions
#' @export
setMethod(
    f = "stripTranscriptVersions",
    signature = signature("Matrix"),
    definition = `stripTranscriptVersions,Matrix`
)



## Updated 2019-07-22.
`stripTranscriptVersions,SummarizedExperiment` <-  # nolint
    `stripTranscriptVersions,matrix`



#' @rdname stripTranscriptVersions
#' @export
setMethod(
    f = "stripTranscriptVersions",
    signature = signature("SummarizedExperiment"),
    definition = `stripTranscriptVersions,SummarizedExperiment`
)



#' @rdname stripTranscriptVersions
#' @export
stripGeneVersions <- function(...) {
    stripTranscriptVersions(...)
}
