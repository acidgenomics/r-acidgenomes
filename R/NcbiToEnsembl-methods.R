#' @name NcbiToEnsembl
#' @inherit AcidGenerics::NcbiToEnsembl description return title
#' @note Updated 2023-11-28.
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
        assert(isOrganism(organism))
        switch(
            EXPR = organism,
            "Homo sapiens" = {
                alert(sprintf(
                    "Matching %d %s against HGNC database.",
                    length(object),
                    ngettext(
                        n = length(object),
                        msg1 = "identifier",
                        msg2 = "identifiers"
                    )
                ))
                hgnc <- Hgnc()
                df <- NcbiToEnsembl(hgnc)
            },
            "Mus musculus" = {
                alert(sprintf(
                    "Matching %d %s against MGI database.",
                    length(object),
                    ngettext(
                        n = length(object),
                        msg1 = "identifier",
                        msg2 = "identifiers"
                    )
                ))
                mgi <- Mgi()
                df <- NcbiToEnsembl(mgi)
            },
            {
                df <- .getEnsemblToNcbiFromOrgDb(
                    keys = as.character(object),
                    organism = organism,
                    return = "NcbiToEnsembl"
                )
                df <- .makeEnsemblToNcbi(
                    object = df,
                    return = "NcbiToEnsembl"
                )
            }
        )
        i <- match(x = object, table = df[[1L]])
        if (anyNA(i)) {
            fail <- object[is.na(i)]
            abort(sprintf(
                "%d match %s: %s.",
                length(fail),
                ngettext(
                    n = length(fail),
                    msg1 = "failure",
                    msg2 = "failures"
                ),
                toInlineString(as.character(fail), n = 10L)
            ))
        }
        out <- df[i, , drop = FALSE]
        rownames(out) <- unname(object)
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
        df <- decode(df)
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


## Updated 2023-11-27.
`NcbiToEnsembl,Mgi` <- # nolint
    `NcbiToEnsembl,Hgnc`


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
    signature = signature(object = "Mgi"),
    definition = `NcbiToEnsembl,Mgi`
)

#' @rdname NcbiToEnsembl
#' @export
setMethod(
    f = "NcbiToEnsembl",
    signature = signature(object = "integer"),
    definition = `NcbiToEnsembl,integer`
)
