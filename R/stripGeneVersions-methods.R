#' @name stripGeneVersions
#' @inherit AcidGenerics::stripGeneVersions
#' @note Updated 2021-08-03.
#'
#' @inheritParams AcidRoxygen::params
#' @param ... Additional arguments.
#'
#' @note This method is strict, and will only strip Ensembl gene identifiers
#'   beginning with `"ENS.*G"`.
#'
#' @seealso
#' - http://www.ensembl.org/info/genome/stable_ids/index.html
#'
#' @examples
#' ## Ensembl.
#' stripGeneVersions("ENSG00000002586.20")
#'
#' ## GENCODE.
#' stripGeneVersions("ENSG00000002586.20_PAR_Y")
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



## Updated 2021-01-27.
`stripGeneVersions,matrix` <-  # nolint
    function(object) {
        assert(hasRownames(object))
        rownames <- rownames(object)
        rownames <- stripGeneVersions(rownames)
        rownames(object) <- rownames
        object
    }



## Updated 2021-01-27.
`stripGeneVersions,Matrix` <-  # nolint
    `stripGeneVersions,matrix`



#' @rdname stripGeneVersions
#' @export
setMethod(
    f = "stripGeneVersions",
    signature = signature("Matrix"),
    definition = `stripGeneVersions,Matrix`
)

#' @rdname stripGeneVersions
#' @export
setMethod(
    f = "stripGeneVersions",
    signature = signature("character"),
    definition = `stripGeneVersions,character`
)

#' @rdname stripGeneVersions
#' @export
setMethod(
    f = "stripGeneVersions",
    signature = signature("matrix"),
    definition = `stripGeneVersions,matrix`
)
