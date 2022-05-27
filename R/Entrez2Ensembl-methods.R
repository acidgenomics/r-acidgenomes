#' @name Entrez2Ensembl
#' @inherit AcidGenerics::Entrez2Ensembl description return title
#' @note Updated 2022-05-27.
#'
#' @inheritParams Ensembl2Entrez
#' @param ... Additional arguments.
#'
#' @examples
#' ## integer ====
#' genes <- c(1L, 2L)
#' organism <- "Homo sapiens"
#' x <- Entrez2Ensembl(object = genes, organism = organism)
#' print(x)
NULL



## Updated 2022-05-27.
`Entrez2Ensembl,integer` <- # nolint
    function(object, organism, format) {
        df <- .getEnsembl2EntrezFromOrgDb(
            keys = as.character(object),
            keytype = "ENTREZID",
            columns = "ENSEMBL",
            organism = organism
        )
        out <- .makeEnsembl2Entrez(
            object = df,
            format = match.arg(format),
            return = "Entrez2Ensembl"
        )
        if (identical(format, "1:1")) {
            idx <- match(x = object, table = out[[1L]])
            out <- out[idx, ]
        }
        out
    }

formals(`Entrez2Ensembl,integer`)[["format"]] <- # nolint
    formals(.makeEnsembl2Entrez)[["format"]]



#' @rdname Entrez2Ensembl
#' @export
setMethod(
    f = "Entrez2Ensembl",
    signature = signature(object = "integer"),
    definition = `Entrez2Ensembl,integer`
)
