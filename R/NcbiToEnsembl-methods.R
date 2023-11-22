#' @name NcbiToEnsembl
#' @inherit AcidGenerics::NcbiToEnsembl description return title
#' @note Updated 2023-11-22.
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
        format <- match.arg(format)
        df <- .getEnsemblToNcbiFromOrgDb(
            keys = as.character(object),
            keytype = "ENTREZID",
            columns = "ENSEMBL",
            organism = organism,
            strict = strict
        )
        assert(identical(object, unique(df[[1L]])))
        out <- .makeEnsemblToNcbi(
            object = df,
            format = format,
            strict = strict,
            return = "NcbiToEnsembl"
        )
        assert(areSetEqual(object, unique(out[[1L]])))
        if (identical(format, "1:1")) {
            i <- match(x = object, table = out[[1L]])
            out <- out[i, , drop = FALSE]
        }
        out
    }

formals(`NcbiToEnsembl,integer`)[["format"]] <- # nolint
    formals(.makeEnsemblToNcbi)[["format"]]



## Updated 2023-11-22.
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
        out <- .makeEnsemblToNcbi(
            object = df,
            format = "1:1",
            strict = TRUE,
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
