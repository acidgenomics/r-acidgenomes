#' @name EnsemblToNcbi
#' @inherit AcidGenerics::EnsemblToNcbi description return title
#' @note Updated 2023-11-21.
#'
#' @inheritParams AcidRoxygen::params
#' @param ... Additional arguments.
#'
#' @param format `character(1)`.
#' Formatting method to apply:
#'
#' - `"1:1"`: *Recommended.* Return with 1:1 mappings. For Ensembl genes that
#' don't map 1:1 with NCBI, pick the oldest NCBI identifier.
#' - `"long"`: Return `1:many` in long format.
#'
#' @param strict `logical(1)`.
#' Error on any mismatches, otherwise return `NA`.
#'
#' @examples
#' ## character ====
#' x <- EnsemblToNcbi(
#'     object = c("ENSG00000000003", "ENSG00000000005"),
#'     organism = "Homo sapiens"
#' )
#' print(x)
NULL



#' Make an `EnsemblToNcbi` (or `NcbiToEnsembl`) object
#'
#' @note Updated 2023-11-21.
#' @noRd
.makeEnsemblToNcbi <-
    function(object,
             format = c("1:1", "long"),
             return = c("EnsemblToNcbi", "NcbiToEnsembl"),
             strict = TRUE) {
        assert(isFlag(strict))
        format <- match.arg(format)
        return <- match.arg(return)
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
        assert(
            is(object, "DFrame"),
            hasRows(object),
            isSubset(cols, colnames(object)),
            !anyNA(object[[fromCol]])
        )
        df <- object[, cols, drop = FALSE]
        df <- decode(df)
        switch(
            EXPR = format,
            "1:1" = {
                if (isAny(df[[2L]], c("List", "list"))) {
                    if (nrow(df) >= 100L) {
                        alert("Mapping nested values 1:1.")
                    }
                    x <- mclapply(
                        X = df[, 2L],
                        FUN = function(x) {
                            x <- sort(
                                x = x,
                                decreasing = FALSE,
                                na.last = TRUE
                            )
                            x <- x[[1L]]
                            x
                        }
                    )
                    x <- unlist(x, recursive = FALSE, use.names = FALSE)
                    df[[2L]] <- x
                } else if (hasDuplicates(df[[1L]])) {
                    i <- order(df, decreasing = FALSE, na.last = TRUE)
                    df <- df[i, , drop = FALSE]
                    i <- !duplicated(df[[1L]])
                    df <- df[i, , drop = FALSE]
                }
                if (!hasRownames(object)) {
                    rownames(df) <- df[[1L]]
                }
                assert(hasNoDuplicates(df[[1L]]))
            },
            "long" = {
                rownames(df) <- NULL
                df <- expand(df)
                i <- order(df, decreasing = FALSE, na.last = TRUE)
                df <- df[i, , drop = FALSE]
            }
        )
        if (isTRUE(strict)) {
            df <- df[complete.cases(df), , drop = FALSE]
        }
        metadata(df) <- append(
            x = metadata(object),
            values = list(
                "date" = Sys.Date(),
                "format" = format,
                "packageVersion" = .pkgVersion,
                "strict" = strict
            )
        )
        new(Class = return, df)
    }



## Updated 2023-11-21.
`EnsemblToNcbi,character` <- # nolint
    function(object,
             organism = NULL,
             format,
             strict = TRUE) {
        assert(
            isCharacter(object),
            hasNoDuplicates(object)
        )
        format <- match.arg(format)
        if (is.null(organism)) {
            organism <- detectOrganism(object)
        }
        if (allAreMatchingFixed(x = object, pattern = ".")) {
            object <- stripGeneVersions(object)
        }
        df <- .getEnsemblToNcbiFromOrgDb(
            keys = object,
            keytype = "ENSEMBL",
            columns = "ENTREZID",
            organism = organism,
            strict = strict
        )
        assert(identical(object, unique(df[[1L]])))
        out <- .makeEnsemblToNcbi(
            object = df,
            format = format,
            return = "EnsemblToNcbi",
            strict = strict
        )
        assert(areSetEqual(object, unique(out[[1L]])))
        out
    }

formals(`EnsemblToNcbi,character`)[["format"]] <- # nolint
    formals(.makeEnsemblToNcbi)[["format"]]



## Updated 2023-11-21.
`EnsemblToNcbi,EnsemblGenes` <- # nolint
    function(object, format) {
        assert(validObject(object))
        df <- mcols(object)
        colnames(df)[colnames(df) == "geneId"] <- "ensemblGeneId"
        metadata(df) <- metadata(object)
        out <- .makeEnsemblToNcbi(
            object = df,
            format = match.arg(format),
            strict = TRUE
        )
        out
    }

formals(`EnsemblToNcbi,EnsemblGenes`)[["format"]] <- # nolint
    formals(.makeEnsemblToNcbi)[["format"]]



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
    signature = signature(object = "character"),
    definition = `EnsemblToNcbi,character`
)
