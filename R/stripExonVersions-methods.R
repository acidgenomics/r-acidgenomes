#' @name stripExonVersions
#' @inherit AcidGenerics::stripExonVersions
#' @note Updated 2024-01-02.
#'
#' @inheritParams AcidRoxygen::params
#' @param ... Additional arguments.
#'
#' @note This method is strict, and will only strip Ensembl exon identifiers
#' beginning with `"ENS.*G"`.
#'
#' @seealso
#' - https://www.ensembl.org/info/genome/stable_ids/index.html
#'
#' @examples
#' ## Ensembl.
#' stripExonVersions("ENSG00000002586.20")
#'
#' ## GENCODE.
#' stripExonVersions("ENSG00000116288.13_10")
#' stripExonVersions("ENSG00000002586.20_PAR_Y")
NULL



## Updated 2023-11-28.
`stripExonVersions,character` <- # nolint
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
`stripExonVersions,integer` <- # nolint
    function(object) {
        object
    }



## Updated 2019-07-22.
`stripExonVersions,matrix` <- # nolint
    function(object) {
        assert(hasRownames(object))
        rownames <- rownames(object)
        rownames <- stripExonVersions(rownames)
        rownames(object) <- rownames
        object
    }



## Updated 2020-01-30.
`stripExonVersions,Matrix` <- # nolint
    `stripExonVersions,matrix`



#' @rdname stripExonVersions
#' @export
setMethod(
    f = "stripExonVersions",
    signature = signature(object = "Matrix"),
    definition = `stripExonVersions,Matrix`
)

#' @rdname stripExonVersions
#' @export
setMethod(
    f = "stripExonVersions",
    signature = signature(object = "character"),
    definition = `stripExonVersions,character`
)

#' @rdname stripExonVersions
#' @export
setMethod(
    f = "stripExonVersions",
    signature = signature(object = "integer"),
    definition = `stripExonVersions,integer`
)

#' @rdname stripExonVersions
#' @export
setMethod(
    f = "stripExonVersions",
    signature = signature(object = "matrix"),
    definition = `stripExonVersions,matrix`
)
