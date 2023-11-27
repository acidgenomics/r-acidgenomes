#' @name NcbiToEnsembl
#' @inherit AcidGenerics::NcbiToEnsembl description return title
#' @note Updated 2023-11-27.
#'
#' @inheritParams EnsemblToNcbi
#' @param ... Additional arguments.
#'
#' @examples
#' ## integer ====
#' x <- NcbiToEnsembl(object = c(2L, 1L), organism = "Homo sapiens")
#' print(x)
NULL



## Updated 2023-11-27.
`NcbiToEnsembl,integer` <- # nolint
    function(object, organism) {
        df <- .getEnsemblToNcbiFromOrgDb(
            keys = as.character(object),
            organism = organism,
            return = "NcbiToEnsembl"
        )
        out <- .makeEnsemblToNcbi(
            object = df,
            return = "NcbiToEnsembl"
        )
        i <- match(x = object, table = out[[1L]])
        assert(!anyNA(i))
        out <- out[i, , drop = FALSE]
        out
    }



## Updated 2023-11-27.
`NcbiToEnsembl,Hgnc` <- # nolint
    function(object) {
        j <- c("ncbiGeneId", "ensemblGeneId")
        assert(
            validObject(object),
            isSubset(j, colnames(object))
        )
        df <- as(object, "DFrame")
        df <- df[, j, drop = FALSE]
        i <- complete.cases(df)
        df <- df[i, , drop = FALSE]
        i <- order(df)
        df <- df[i, , drop = FALSE]
        i <- !duplicated(df[[1L]]) & !duplicated(df[[2L]])
        df <- df[i, , drop = FALSE]
        out <- .makeEnsemblToNcbi(
            object = df,
            return = "NcbiToEnsembl"
        )
        out
    }



#' @rdname NcbiToEnsembl
#' @export
setMethod(
    f = "NcbiToEnsembl",
    signature = signature(object = "Hgnc"),
    definition = `NcbiToEnsembl,Hgnc`
)

#' @rdname NcbiToEnsembl
#' @export
setMethod(
    f = "NcbiToEnsembl",
    signature = signature(object = "integer"),
    definition = `NcbiToEnsembl,integer`
)
