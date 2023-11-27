#' @name EnsemblToNcbi
#' @inherit AcidGenerics::EnsemblToNcbi description return title
#' @note Updated 2023-11-27.
#'
#' @inheritParams AcidRoxygen::params
#' @param ... Additional arguments.
#'
#' @param strict `logical(1)`.
#' Error on any mismatches, otherwise return `NA`.
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
    function(object, strict, return) {
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
            isFlag(strict),
            isOrganism(organism),
            isSubset(cols, colnames(object)),
            !anyNA(object[[fromCol]])
        )
        df <- object[, cols, drop = FALSE]
        df <- decode(df)
        df <- expand(df)
        i <- order(df, decreasing = FALSE, na.last = TRUE)
        df <- df[i, , drop = FALSE]
        if (organism == "Homo sapiens" && anyDuplicated(df[[2L]])) {
            alert("Resolving ambiguous duplicates with HGNC annotations.")
            hgnc <- Hgnc()
            df2 <- EnsemblToNcbi(hgnc)
            df2 <- df2[, cols, drop = FALSE]
            i <- df2[[1L]] %in% df[[1L]]
            df2 <- df2[i, , drop = FALSE]
            i <- order(df2, decreasing = FALSE, na.last = TRUE)
            df2 <- df2[i, , drop = FALSE]
            df <- rbind(df2, df)
        }
        i <- !duplicated(df[[1L]]) & !duplicated(df[[2L]])
        df <- df[i, , drop = FALSE]
        i <- order(df, decreasing = FALSE, na.last = TRUE)
        df <- df[i, , drop = FALSE]
        assert(
            hasNoDuplicates(df[[1L]]),
            hasNoDuplicates(df[[2L]])
        )
        ## FIXME Rework NA filling in strict mode.
        if (isTRUE(strict)) {
            i <- complete.cases(df)
            df <- df[i, , drop = FALSE]
        }
        metadata(df) <- append(
            x = metadata(object),
            values = list(
                "date" = Sys.Date(),
                "format" = "1:1",
                "organism" = organism,
                "packageVersion" = .pkgVersion,
                "strict" = strict
            )
        )
        new(Class = return, df)
    }



## Updated 2023-11-27.
`EnsemblToNcbi,character` <- # nolint
    function(object,
             organism = NULL,
             strict = TRUE) {
        assert(
            isCharacter(object),
            hasNoDuplicates(object)
        )
        keys <- unname(object)
        if (is.null(organism)) {
            organism <- detectOrganism(object)
        }
        if (allAreMatchingFixed(x = keys, pattern = ".")) {
            keys <- stripGeneVersions(keys)
        }
        df <- .getEnsemblToNcbiFromOrgDb(
            keys = keys,
            organism = organism,
            return = "EnsemblToNcbi"
        )
        out <- .makeEnsemblToNcbi(
            object = df,
            strict = strict,
            return = "EnsemblToNcbi"
        )
        i <- match(x = keys, table = out[[1L]])
        assert(!anyNA(i))
        out <- out[i, , drop = FALSE]
        rownames(out) <- unname(object)
        out
    }



## Updated 2023-11-27.
`EnsemblToNcbi,EnsemblGenes` <- # nolint
    function(object) {
        assert(validObject(object))
        df <- mcols(object)
        colnames(df)[colnames(df) == "geneId"] <- "ensemblGeneId"
        metadata(df) <- metadata(object)
        out <- .makeEnsemblToNcbi(
            object = df,
            strict = TRUE,
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
        i <- complete.cases(df)
        df <- df[i, , drop = FALSE]
        i <- order(df)
        df <- df[i, , drop = FALSE]
        i <- !duplicated(df[[1L]]) & !duplicated(df[[2L]])
        df <- df[i, , drop = FALSE]
        out <- .makeEnsemblToNcbi(
            object = df,
            strict = TRUE,
            return = "EnsemblToNcbi"
        )
        out
    }



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
    signature = signature(object = "character"),
    definition = `EnsemblToNcbi,character`
)
