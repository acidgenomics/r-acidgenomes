#' @name EnsemblToNcbi
#' @inherit AcidGenerics::EnsemblToNcbi description return title
#' @note Updated 2023-11-27.
#'
#' @inheritParams AcidRoxygen::params
#' @param ... Additional arguments.
#'
#' @examples
#' ## character ====
#' x <- EnsemblToNcbi(
#'     object = c("ENSG00000000005.6", "ENSG00000000003.16"),
#'     organism = "Homo sapiens"
#' )
#' print(x)
NULL



#' Make an `EnsemblToNcbi` (or `NcbiToEnsembl`) object
#'
#' @note Updated 2023-11-27.
#' @noRd
.makeEnsemblToNcbi <-
    function(object, return) {
        return <- match.arg(
            arg = return,
            choices = c("EnsemblToNcbi", "NcbiToEnsembl")
        )
        switch(
            EXPR = return,
            "EnsemblToNcbi" = {
                fromCol <- "ensemblGeneId"
                toCol <- "ncbiGeneId"
            },
            "NcbiToEnsembl" = {
                fromCol <- "ncbiGeneId"
                toCol <- "ensemblGeneId"
            }
        )
        cols <- c(fromCol, toCol)
        organism <- metadata(object)[["organism"]]
        assert(
            is(object, "DFrame"),
            hasRows(object),
            isOrganism(organism),
            isSubset(cols, colnames(object)),
            !anyNA(object[[fromCol]])
        )
        map <- object[, cols, drop = FALSE]
        map <- decode(map)
        map <- expand(map)
        i <- complete.cases(map)
        map <- map[i, , drop = FALSE]
        i <- order(map)
        map <- map[i, , drop = FALSE]
        i <- !duplicated(map[[1L]]) & !duplicated(map[[2L]])
        map <- map[i, , drop = FALSE]
        i <- order(map)
        map <- map[i, , drop = FALSE]
        assert(
            all(complete.cases(map)),
            hasNoDuplicates(map[[1L]]),
            hasNoDuplicates(map[[2L]])
        )
        metadata(map) <- append(
            x = metadata(object),
            values = list(
                "date" = Sys.Date(),
                "organism" = organism,
                "packageVersion" = .pkgVersion
            )
        )
        new(Class = return, map)
    }



## Updated 2023-11-27.
`EnsemblToNcbi,character` <- # nolint
    function(object, organism = NULL) {
        if (is.null(organism)) {
            organism <- detectOrganism(object)
        }
        keys <- unname(object)
        if (allAreMatchingFixed(x = keys, pattern = ".")) {
            keys <- stripGeneVersions(keys)
        }
        switch(
            EXPR = organism,
            "Homo sapiens" = {
                hgnc <- Hgnc()
                df <- EnsemblToNcbi(hgnc)
            },
            "Mus musculus" = {
                mgi <- Mgi()
                df <- EnsemblToNcbi(mgi)
            },
            {
                df <- .getEnsemblToNcbiFromOrgDb(
                    keys = keys,
                    organism = organism,
                    return = "EnsemblToNcbi"
                )
                df <- .makeEnsemblToNcbi(
                    object = df,
                    return = "EnsemblToNcbi"
                )
            }
        )
        ## FIXME Use Hgnc or Mgi methods instead for Homo sapiens and Mus musculus.
        i <- match(x = keys, table = df[[1L]])
        assert(!anyNA(i))
        out <- df[i, , drop = FALSE]
        rownames(out) <- unname(object)
        out
    }



## FIXME Use Hgnc to pick curated matches for Homo sapiens.
## FIXME Use Mgi to pick curated matches for Mus musculus.

## Updated 2023-11-27.
`EnsemblToNcbi,EnsemblGenes` <- # nolint
    function(object) {
        assert(validObject(object))
        df <- mcols(object)
        colnames(df)[colnames(df) == "geneId"] <- "ensemblGeneId"
        metadata(df) <- metadata(object)
        out <- .makeEnsemblToNcbi(
            object = df,
            return = "EnsemblToNcbi"
        )
        out
    }



## Updated 2023-11-27.
`EnsemblToNcbi,Hgnc` <- # nolint
    function(object) {
        j <- c("ensemblGeneId", "ncbiGeneId")
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
            return = "EnsemblToNcbi"
        )
        out
    }



## Updated 2023-11-27.
`EnsemblToNcbi,Mgi` <- # nolint
    `EnsemblToNcbi,Hgnc`



## Updated 2023-03-01.
`EnsemblToNcbi,GencodeGenes` <- # nolint
    `EnsemblToNcbi,EnsemblGenes`



#' @rdname EnsemblToNcbi
#' @export
setMethod(
    f = "EnsemblToNcbi",
    signature = signature(object = "EnsemblGenes"),
    definition = `EnsemblToNcbi,EnsemblGenes`
)

#' @rdname EnsemblToNcbi
#' @export
setMethod(
    f = "EnsemblToNcbi",
    signature = signature(object = "GencodeGenes"),
    definition = `EnsemblToNcbi,GencodeGenes`
)

#' @rdname EnsemblToNcbi
#' @export
setMethod(
    f = "EnsemblToNcbi",
    signature = signature(object = "Hgnc"),
    definition = `EnsemblToNcbi,Hgnc`
)

#' @rdname EnsemblToNcbi
#' @export
setMethod(
    f = "EnsemblToNcbi",
    signature = signature(object = "Mgi"),
    definition = `EnsemblToNcbi,Mgi`
)

#' @rdname EnsemblToNcbi
#' @export
setMethod(
    f = "EnsemblToNcbi",
    signature = signature(object = "character"),
    definition = `EnsemblToNcbi,character`
)
