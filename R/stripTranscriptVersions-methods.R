#' @name stripTranscriptVersions
#' @inherit AcidGenerics::stripTranscriptVersions
#' @note Updated 2021-08-03.
#'
#' @inheritParams AcidRoxygen::params
#' @param ... Additional arguments.
#'
#' @note This method is strict, and will only strip Ensembl transcript
#'   identifiers beginning with `"ENS.*T"`.
#'
#' @seealso
#' - http://www.ensembl.org/info/genome/stable_ids/index.html
#'
#' @examples
#' ## Ensembl.
#' stripTranscriptVersions("ENST00000431238.7")
#'
#' ## GENCODE.
#' stripTranscriptVersions("ENST00000431238.7_PAR_Y")
#'
#' ## WormBase (no modification).
#' stripTranscriptVersions("cTel79B.1")
NULL



## Updated 2021-01-27.
`stripTranscriptVersions,character` <-  # nolint
    function(object) {
        assert(isCharacter(object))
        pattern <- "^(ENS.*T[0-9]{11})(\\.[0-9]+)(_.+)?$"
        if (!any(grepl(pattern = pattern, x = object))) {
            alertWarning("No transcript versions to modify.")
            return(object)
        }
        out <- gsub(
            pattern = pattern,
            replacement = "\\1\\3",
            x = object
        )
        out
    }



## Updated 2022-03-09.
`stripTranscriptVersions,integer` <-  # nolint
    function(object) {
        return(object)
    }



## Updated 2019-07-22.
`stripTranscriptVersions,matrix` <-  # nolint
    function(object) {
        assert(hasRownames(object))
        rownames <- rownames(object)
        rownames <- stripTranscriptVersions(rownames)
        rownames(object) <- rownames
        object
    }



## Updated 2020-01-30.
`stripTranscriptVersions,Matrix` <-  # nolint
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
