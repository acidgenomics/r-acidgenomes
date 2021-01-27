## FIXME CONSIDER ADDING SUPPORT FOR GRANGES.
##       IF ADDED, ALSO MODIFY MCOLS IN ROWDATA FOR SE.



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
#' @seealso
#' - http://www.ensembl.org/info/genome/stable_ids/index.html
#'
#' @examples
#' ## Ensembl.
#' stripGeneVersions("ENSG00000002586.20")
#' stripTranscriptVersions("ENST00000431238.7")
#'
#' ## GENCODE.
#' stripGeneVersions("ENSG00000002586.20_PAR_Y")
#' stripTranscriptVersions("ENST00000431238.7_PAR_Y")
#'
#' ## WormBase.
#' stripTranscriptVersions("cTel79B.1")
NULL




## Updated 2021-01-27.
`stripGeneVersions,character` <-  # nolint
    function(object) {
        assert(isCharacter(object))
        pattern <- "^(ENS.*G[0-9]{11})(\\.[0-9]+)(_.+)?$"
        if (!any(grepl(pattern = pattern, x = object))) {
            alertWarning("No gene versions to modify.")
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



#' @rdname stripVersions
#' @export
setMethod(
    f = "stripTranscriptVersions",
    signature = signature("character"),
    definition = `stripTranscriptVersions,character`
)



## Updated 2021-01-27.
`stripGeneVersions,matrix` <-  # nolint
    function(object) {
        assert(hasRownames(object))
        rownames <- rownames(object)
        rownames <- stripGeneVersions(rownames)
        rownames(object) <- rownames
        object
    }



#' @rdname stripVersions
#' @export
setMethod(
    f = "stripGeneVersions",
    signature = signature("matrix"),
    definition = `stripGeneVersions,matrix`
)



## Updated 2019-07-22.
`stripTranscriptVersions,matrix` <-  # nolint
    function(object) {
        assert(hasRownames(object))
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



## Updated 2021-01-27.
`stripGeneVersions,Matrix` <-  # nolint
    `stripGeneVersions,matrix`



#' @rdname stripVersions
#' @export
setMethod(
    f = "stripGeneVersions",
    signature = signature("Matrix"),
    definition = `stripGeneVersions,Matrix`
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
`stripGeneVersions,SummarizedExperiment` <-  # nolint
    `stripGeneVersions,matrix`



#' @rdname stripVersions
#' @export
setMethod(
    f = "stripGeneVersions",
    signature = signature("SummarizedExperiment"),
    definition = `stripGeneVersions,SummarizedExperiment`
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
