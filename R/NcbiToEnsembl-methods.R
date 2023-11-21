#' @name NcbiToEnsembl
#' @inherit AcidGenerics::NcbiToEnsembl description return title
#' @note Updated 2023-09-16.
#'
#' @inheritParams EnsemblToNcbi
#' @param ... Additional arguments.
#'
#' @examples
#' ## integer ====
#' x <- NcbiToEnsembl(object = c(1L, 2L), organism = "Homo sapiens")
#' print(x)
NULL



## Updated 2023-11-21.
`NcbiToEnsembl,integer` <- # nolint
    function(object, organism, format, strict = TRUE) {
        df <- .getEnsemblToNcbiFromOrgDb(
            keys = as.character(object),
            keytype = "ENTREZID",
            columns = "ENSEMBL",
            organism = organism,
            strict = strict
        )
        out <- .makeEnsemblToNcbi(
            object = df,
            format = match.arg(format),
            return = "NcbiToEnsembl",
            strict = strict
        )
        if (identical(format, "1:1")) {
            idx <- match(x = object, table = out[[1L]])
            out <- out[idx, , drop = FALSE]
        }
        out
    }

formals(`NcbiToEnsembl,integer`)[["format"]] <- # nolint
    formals(.makeEnsemblToNcbi)[["format"]]



#' @rdname NcbiToEnsembl
#' @export
setMethod(
    f = "NcbiToEnsembl",
    signature = signature(object = "integer"),
    definition = `NcbiToEnsembl,integer`
)
