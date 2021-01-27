#' @name stripVersions
#' @inherit AcidGenerics::stripTranscriptVersions
#' @note Updated 2021-01-27.
#'
#' @inheritParams AcidRoxygen::params
#' @param ... Additional arguments.
#'
#' @note This method is strict, and will only strip Ensembl transcript IDs
#'   beginning with "ENS".
#'
#' @examples
#' ## Ensembl.
#' stripTranscriptVersions("ENSMUST00000000001.1")
#'
#' ## GENCODE.
#' stripGeneVersions("ENSG00000002586.20_PAR_Y")
#'
#' ## WormBase.
#' stripTranscriptVersions("cTel79B.1")
NULL




## Updated 2021-01-27.
`stripGeneVersions,character` <-  # nolint
    function(object) {
        assert(isCharacter(object))
        pattern <- "^(ENS.*G[0-9]{11})(\\.[0-9]+)(_.+)?$"
        if (!any(isMatchingRegex(pattern = pattern, x = object))) {
            alertWarning("No gene versions detected.")
            return(object)
        }
        out <- gsub(
            pattern = pattern,
            replacement = "\\1\\3",
            x = object
        )
        out
    }



#' @rdname stripVersions
#' @export
setMethod(
    f = "stripGeneVersions",
    signature = signature("character"),
    definition = `stripGeneVersions,character`
)




## FIXME Make this transcript specific, see gene code above.


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



#' @rdname stripVersions
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



#' @rdname stripVersions
#' @export
setMethod(
    f = "stripTranscriptVersions",
    signature = signature("matrix"),
    definition = `stripTranscriptVersions,matrix`
)



## Updated 2020-01-30.
`stripTranscriptVersions,Matrix` <-  # nolint
    `stripTranscriptVersions,matrix`



#' @rdname stripVersions
#' @export
setMethod(
    f = "stripTranscriptVersions",
    signature = signature("Matrix"),
    definition = `stripTranscriptVersions,Matrix`
)



## Updated 2019-07-22.
`stripTranscriptVersions,SummarizedExperiment` <-  # nolint
    `stripTranscriptVersions,matrix`



#' @rdname stripVersions
#' @export
setMethod(
    f = "stripTranscriptVersions",
    signature = signature("SummarizedExperiment"),
    definition = `stripTranscriptVersions,SummarizedExperiment`
)
