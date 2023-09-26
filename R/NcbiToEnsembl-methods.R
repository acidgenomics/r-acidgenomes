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



## Updated 2022-05-27.
`NcbiToEnsembl,integer` <- # nolint
    function(object, organism, format) {
        df <- .getEnsemblToNcbiFromOrgDb(
            keys = as.character(object),
            keytype = "ENTREZID",
            columns = "ENSEMBL",
            organism = organism
        )
        out <- .makeEnsemblToNcbi(
            object = df,
            format = match.arg(format),
            return = "NcbiToEnsembl"
        )
        if (identical(format, "1:1")) {
            idx <- match(x = object, table = out[[1L]])
            out <- out[idx, ]
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
