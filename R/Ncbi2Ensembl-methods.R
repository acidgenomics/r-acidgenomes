#' @name Ncbi2Ensembl
#' @inherit AcidGenerics::Ncbi2Ensembl description return title
#' @note Updated 2023-03-01.
#'
#' @inheritParams Ensembl2Ncbi
#' @param ... Additional arguments.
#'
#' @examples
#' ## integer ====
#' genes <- c(1L, 2L)
#' organism <- "Homo sapiens"
#' x <- Ncbi2Ensembl(object = genes, organism = organism)
#' print(x)
NULL



## Updated 2022-05-27.
`Ncbi2Ensembl,integer` <- # nolint
    function(object, organism, format) {
        df <- .getEnsembl2NcbiFromOrgDb(
            keys = as.character(object),
            keytype = "ENTREZID",
            columns = "ENSEMBL",
            organism = organism
        )
        out <- .makeEnsembl2Ncbi(
            object = df,
            format = match.arg(format),
            return = "Ncbi2Ensembl"
        )
        if (identical(format, "1:1")) {
            idx <- match(x = object, table = out[[1L]])
            out <- out[idx, ]
        }
        out
    }

formals(`Ncbi2Ensembl,integer`)[["format"]] <- # nolint
    formals(.makeEnsembl2Ncbi)[["format"]]



#' @rdname Ncbi2Ensembl
#' @export
setMethod(
    f = "Ncbi2Ensembl",
    signature = signature(object = "integer"),
    definition = `Ncbi2Ensembl,integer`
)
