#' @name stripTranscriptVersions
#' @inherit AcidGenerics::stripTranscriptVersions
#' @note Updated 2023-11-28.
#'
#' @inheritParams AcidRoxygen::params
#' @param ... Additional arguments.
#'
#' @note This method is strict, and will only strip Ensembl transcript
#' identifiers beginning with `"ENS.*T"`.
#'
#' @seealso
#' - https://www.ensembl.org/info/genome/stable_ids/index.html
#'
#' @examples
#' ## Ensembl.
#' stripTranscriptVersions("ENST00000431238.7")
#'
#' ## GENCODE.
#' stripTranscriptVersions("ENST00000493373.5_7")
#' stripTranscriptVersions("ENST00000431238.7_PAR_Y")
#'
#' ## WormBase (no modification).
#' stripTranscriptVersions("cTel79B.1")
NULL



## Updated 2023-11-28.
`stripTranscriptVersions,character` <- # nolint
    function(object) {
        assert(isCharacter(object))
        pattern <- "^(ENS.*T[0-9]{11})(\\.[0-9]+)(_[0-9]+)?(_PAR_Y)?$"
        if (!any(grepl(pattern = pattern, x = object))) {
            alertWarning("No transcript versions to modify.")
            return(object)
        }
        out <- gsub(
            pattern = pattern,
            replacement = "\\1\\4",
            x = object
        )
        out
    }



## Updated 2022-03-09.
`stripTranscriptVersions,integer` <- # nolint
    function(object) {
        object
    }



## Updated 2019-07-22.
`stripTranscriptVersions,matrix` <- # nolint
    function(object) {
        assert(hasRownames(object))
        rownames <- rownames(object)
        rownames <- stripTranscriptVersions(rownames)
        rownames(object) <- rownames
        object
    }



## Updated 2020-01-30.
`stripTranscriptVersions,Matrix` <- # nolint
    `stripTranscriptVersions,matrix`



#' @rdname stripTranscriptVersions
#' @export
setMethod(
    f = "stripTranscriptVersions",
    signature = signature(object = "Matrix"),
    definition = `stripTranscriptVersions,Matrix`
)

#' @rdname stripTranscriptVersions
#' @export
setMethod(
    f = "stripTranscriptVersions",
    signature = signature(object = "character"),
    definition = `stripTranscriptVersions,character`
)

#' @rdname stripTranscriptVersions
#' @export
setMethod(
    f = "stripTranscriptVersions",
    signature = signature(object = "integer"),
    definition = `stripGeneVersions,integer`
)

#' @rdname stripTranscriptVersions
#' @export
setMethod(
    f = "stripTranscriptVersions",
    signature = signature(object = "matrix"),
    definition = `stripTranscriptVersions,matrix`
)
