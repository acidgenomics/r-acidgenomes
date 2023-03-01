#' @name Ensembl2Ncbi
#' @inherit AcidGenerics::Ensembl2Ncbi description return title
#' @note Updated 2023-03-01.
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
#' organism <- "Homo sapiens"
#'
#' ## character ====
#' ## Ensembl-to-NCBI.
#' genes <- c("ENSG00000000003", "ENSG00000000005")
#' x <- Ensembl2Ncbi(object = genes, organism = organism)
#' print(x)
NULL



#' Make an `Ensembl2Ncbi` (or `Ncbi2Ensembl`) object
#'
#' @note Updated 2023-03-01.
#' @noRd
.makeEnsembl2Ncbi <-
    function(object,
             format = c("1:1", "long"),
             ## Internal-only args:
             return = c("Ensembl2Ncbi", "Ncbi2Ensembl")) {
        format <- match.arg(format)
        return <- match.arg(return)
        switch(
            EXPR = return,
            "Ensembl2Ncbi" = {
                fromCol <- "ensemblGeneId"
                toCol <- "ncbiGeneId"
            },
            "Ncbi2Ensembl" = {
                fromCol <- "ncbiGeneId"
                toCol <- "ensemblGeneId"
            }
        )
        cols <- c(fromCol, toCol)
        assert(
            is(object, "DataFrame"),
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
        assert(is(df, "DataFrame"))
        df <- df[complete.cases(df), , drop = FALSE]
        metadata(df) <- append(
            x = metadata(object),
            values = list("format" = format)
        )
        new(Class = return, df)
    }



## Updated 2022-05-27.
`Ensembl2Ncbi,character` <- # nolint
    function(object,
             organism = NULL,
             format) {
        if (is.null(organism)) {
            organism <- detectOrganism(object)
        }
        if (allAreMatchingFixed(x = object, pattern = ".")) {
            object <- stripGeneVersions(object)
        }
        df <- .getEnsembl2NcbiFromOrgDb(
            keys = object,
            keytype = "ENSEMBL",
            columns = "ENTREZID",
            organism = organism
        )
        out <- .makeEnsembl2Ncbi(
            object = df,
            format = match.arg(format),
            return = "Ensembl2Ncbi"
        )
        if (identical(format, "1:1")) {
            idx <- match(x = object, table = out[[1L]])
            out <- out[idx, ]
        }
        out
    }

formals(`Ensembl2Ncbi,character`)[["format"]] <- # nolint
    formals(.makeEnsembl2Ncbi)[["format"]]



## FIXME Consider reworking this approach.
## FIXME This should only apply if gene identifier contains Ensembl identifiers.
## FIXME Let's just class this against EnsemblGenes or GencodeGenes instead.

## Updated 2023-03-01.
`Ensembl2Ncbi,GenomicRanges` <- # nolint
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
        .makeEnsembl2Ncbi(
            object = df,
            format = match.arg(format)
        )
    }

formals(`Ensembl2Ncbi,GenomicRanges`)[["format"]] <- # nolint
    formals(.makeEnsembl2Ncbi)[["format"]]



#' @rdname Ensembl2Ncbi
#' @export
setMethod(
    f = "Ensembl2Ncbi",
    signature = signature(object = "GenomicRanges"),
    definition = `Ensembl2Ncbi,GenomicRanges`
)

#' @rdname Ensembl2Ncbi
#' @export
setMethod(
    f = "Ensembl2Ncbi",
    signature = signature(object = "character"),
    definition = `Ensembl2Ncbi,character`
)
