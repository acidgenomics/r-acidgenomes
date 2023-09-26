#' @name EnsemblToNcbi
#' @inherit AcidGenerics::EnsemblToNcbi description return title
#' @note Updated 2023-09-16.
#'
#' @inheritParams AcidRoxygen::params
#' @param ... Additional arguments.
#'
#' @param format `character(1)`.
#' Formatting method to apply:
#'
#' - `"1:1"`: *Recommended.* Return with 1:1 mappings. For Ensembl genes that
#' don't map 1:1 with NCBI, pick the oldest NCBI identifier. Genes that don't
#' map to NCBI will contain `NA` in `ncbiGeneId` column.
#' - `"long"`: Return `1:many` in long format.
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
#' @note Updated 2023-03-01.
#' @noRd
.makeEnsemblToNcbi <-
    function(object,
             format = c("1:1", "long"),
             ## Internal-only args:
             return = c("EnsemblToNcbi", "NcbiToEnsembl")) {
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
            isSubset(cols, colnames(object))
        )
        df <- object[, cols, drop = FALSE]
        switch(
            EXPR = format,
            "1:1" = {
                if (isAny(df[[2L]], c("List", "list"))) {
                    x <- lapply(
                        X = df[, 2L],
                        FUN = function(x) {
                            sort(
                                x = x,
                                decreasing = FALSE,
                                na.last = TRUE
                            )[[1L]]
                        }
                    )
                    x <- unlist(x, recursive = FALSE, use.names = FALSE)
                    df[[2L]] <- x
                } else if (hasDuplicates(df[[1L]])) {
                    df <- df[order(df, decreasing = FALSE, na.last = TRUE), ]
                    spl <- split(x = df, f = df[[1L]])
                    spl <- lapply(
                        X = spl,
                        FUN = function(x) {
                            x[1L, ]
                        }
                    )
                    df <- do.call(what = rbind, args = spl)
                }
                ## This step is useful for character method handling.
                if (!hasRownames(object)) {
                    rownames(df) <- df[[1L]]
                }
                assert(hasNoDuplicates(df[[1L]]))
            },
            "long" = {
                rownames(df) <- NULL
                df <- expand(df)
                df <- df[order(df), ]
            }
        )
        assert(is(df, "DFrame"))
        df <- df[complete.cases(df), , drop = FALSE]
        metadata(df) <- append(
            x = metadata(object),
            values = list(
                "date" = Sys.Date(),
                "format" = format,
                "packageVersion" = .pkgVersion
            )
        )
        new(Class = return, df)
    }



## Updated 2022-05-27.
`EnsemblToNcbi,character` <- # nolint
    function(object,
             organism = NULL,
             format) {
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
            organism = organism
        )
        out <- .makeEnsemblToNcbi(
            object = df,
            format = match.arg(format),
            return = "EnsemblToNcbi"
        )
        if (identical(format, "1:1")) {
            idx <- match(x = object, table = out[[1L]])
            out <- out[idx, ]
        }
        out
    }

formals(`EnsemblToNcbi,character`)[["format"]] <- # nolint
    formals(.makeEnsemblToNcbi)[["format"]]



## Updated 2023-03-01.
`EnsemblToNcbi,EnsemblGenes` <- # nolint
    function(object, format) {
        assert(hasColnames(mcols(object)))
        colnames(mcols(object)) <-
            camelCase(
                object = colnames(mcols(object)),
                strict = TRUE
            )
        assert(
            isSubset(
                x = c("geneId", "ncbiGeneId"),
                y = colnames(mcols(object))
            )
        )
        df <- mcols(object)
        colnames(df)[colnames(df) == "geneId"] <- "ensemblGeneId"
        metadata(df) <- metadata(object)
        out <- .makeEnsemblToNcbi(
            object = df,
            format = match.arg(format)
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
