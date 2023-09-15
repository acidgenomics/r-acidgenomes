#' @name stripGeneVersions
#' @inherit AcidGenerics::stripGeneVersions
#' @note Updated 2022-03-09.
#'
#' @inheritParams AcidRoxygen::params
#' @param ... Additional arguments.
#'
#' @note This method is strict, and will only strip Ensembl gene identifiers
#' beginning with `"ENS.*G"`.
#'
#' @seealso
#' - https://www.ensembl.org/info/genome/stable_ids/index.html
#'
#' @examples
#' ## Ensembl.
#' stripGeneVersions("ENSG00000002586.20")
#'
#' ## GENCODE.
#' stripGeneVersions("ENSG00000002586.20_PAR_Y")
NULL



## Updated 2021-01-27.
`stripGeneVersions,character` <- # nolint
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



## Allow passthrough support of NCBI gene identifiers (e.g. "7157", which
## corresponds to "TP53". This method is used in DepMapAnalysis package.
## Updated 2022-03-09.
`stripGeneVersions,integer` <- # nolint
    function(object) {
        object
    }



## Updated 2021-01-27.
`stripGeneVersions,matrix` <- # nolint
    function(object) {
        assert(hasRownames(object))
        rownames <- rownames(object)
        rownames <- stripGeneVersions(rownames)
        rownames(object) <- rownames
        object
    }



## Updated 2021-01-27.
`stripGeneVersions,Matrix` <- # nolint
    `stripGeneVersions,matrix`



#' @rdname stripGeneVersions
#' @export
setMethod(
    f = "stripGeneVersions",
    signature = signature(object = "Matrix"),
    definition = `stripGeneVersions,Matrix`
)

#' @rdname stripGeneVersions
#' @export
setMethod(
    f = "stripGeneVersions",
    signature = signature(object = "character"),
    definition = `stripGeneVersions,character`
)

#' @rdname stripGeneVersions
#' @export
setMethod(
    f = "stripGeneVersions",
    signature = signature(object = "integer"),
    definition = `stripGeneVersions,integer`
)

#' @rdname stripGeneVersions
#' @export
setMethod(
    f = "stripGeneVersions",
    signature = signature(object = "matrix"),
    definition = `stripGeneVersions,matrix`
)
