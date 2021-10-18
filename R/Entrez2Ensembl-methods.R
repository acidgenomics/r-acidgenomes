#' @name Entrez2Ensembl
#' @inherit AcidGenerics::Entrez2Ensembl description return title
#' @note Updated 2021-08-18.
#'
#' @inheritParams AcidRoxygen::params
#'
#' @examples
#' organism <- "Homo sapiens"
#'
#' ## integer ====
#' ## Entrez-to-Ensembl.
#' genes <- c(1L, 2L)
#' x <- Entrez2Ensembl(object = genes, organism = organism)
#' print(x)
NULL



## Updated 2021-01-18.
`Entrez2Ensembl,integer` <-  # nolint
    function(object, organism, format) {
        df <- .getEnsembl2EntrezFromOrgDb(
            keys = as.character(object),
            keytype = "ENTREZID",
            columns = "ENSEMBL",
            organism = organism
        )
        .makeEnsembl2Entrez(
            object = df,
            format = match.arg(format),
            return = "Entrez2Ensembl"
        )
    }

formals(`Entrez2Ensembl,integer`)[["format"]] <-
    formals(.makeEnsembl2Entrez)[["format"]]



#' @rdname Ensembl2Entrez
#' @export
setMethod(
    f = "Entrez2Ensembl",
    signature = signature(object = "integer"),
    definition = `Entrez2Ensembl,integer`
)
