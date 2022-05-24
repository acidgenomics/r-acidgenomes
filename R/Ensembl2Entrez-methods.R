#' @name Ensembl2Entrez
#' @inherit AcidGenerics::Ensembl2Entrez description return title
#' @note Updated 2022-05-04.
#'
#' @inheritParams AcidRoxygen::params
#' @param ... Additional arguments.
#'
#' @param format `character(1)`.
#' Formatting method to apply:
#'
#' - `"1:1"`: *Recommended.* Return with 1:1 mappings. For Ensembl genes that
#' don't map 1:1 with Entrez, pick the oldest Entrez identifier. Genes that
#' don't map to Entrez will contain `NA` in `entrezId` column.
#' - `"long"`: Return `1:many` in long format.
#'
#' @examples
#' organism <- "Homo sapiens"
#'
#' ## character ====
#' ## Ensembl-to-Entrez.
#' genes <- c("ENSG00000000003", "ENSG00000000005")
#' x <- Ensembl2Entrez(object = genes, organism = organism)
#' print(x)
NULL



#' Make an Ensembl2Entrez (or Entrez2Ensembl) object
#'
#' @note Updated 2021-02-24.
#' @noRd
.makeEnsembl2Entrez <-
    function(object,
             format = c("1:1", "long"),
             ## Internal-only args:
             return = c("Ensembl2Entrez", "Entrez2Ensembl")) {
        format <- match.arg(format)
        return <- match.arg(return)
        cols <- switch(
            EXPR = return,
            "Ensembl2Entrez" = c("ensemblId", "entrezId"),
            "Entrez2Ensembl" = c("entrezId", "ensemblId")
        )
        assert(
            is(object, "DataFrame"),
            hasRows(object),
            isSubset(cols, colnames(object))
        )
        df <- object[, cols, drop = FALSE]
        rownames(df) <- NULL
        df <- decode(expand(df))
        assert(
            is.character(df[["ensemblId"]]),
            is.integer(df[["entrezId"]])
        )
        if (identical(format, "1:1")) {
            x <- split(x = df, f = df[[1L]])
            assert(is(x, "SplitDataFrameList"))
            x <- lapply(
                X = x[, 2L],
                FUN = function(x) {
                    sort(x = x, decreasing = FALSE, na.last = TRUE)[[1L]]
                }
            )
            x <- unlist(x, recursive = FALSE, use.names = TRUE)
            df <- DataFrame("a" = names(x), "b" = unname(x))
            rownames(df) <- df[[1L]]
            colnames(df) <- cols
        }
        df <- df[complete.cases(df), , drop = FALSE]
        assert(hasRows(df))
        df[["entrezId"]] <- as.integer(df[["entrezId"]])
        metadata(df) <- metadata(object)
        metadata(df)[["format"]] <- format
        new(Class = return, df)
    }



## Updated 2021-02-10.
`Ensembl2Entrez,character` <- # nolint
    function(object,
             organism = NULL,
             format) {
        if (is.null(organism)) {
            organism <- detectOrganism(object)
        }
        if (allAreMatchingFixed(x = object, pattern = ".")) {
            object <- stripGeneVersions(object)
        }
        df <- .getEnsembl2EntrezFromOrgDb(
            keys = object,
            keytype = "ENSEMBL",
            columns = "ENTREZID",
            organism = organism
        )
        .makeEnsembl2Entrez(
            object = df,
            format = match.arg(format),
            return = "Ensembl2Entrez"
        )
    }

formals(`Ensembl2Entrez,character`)[["format"]] <- # nolint
    formals(.makeEnsembl2Entrez)[["format"]]



## Updated 2021-02-01.
`Ensembl2Entrez,GenomicRanges` <- # nolint
    function(object, format) {
        assert(hasColnames(mcols(object)))
        colnames(mcols(object)) <-
            camelCase(
                object = colnames(mcols(object)),
                strict = TRUE
            )
        assert(
            isSubset(
                x = c("geneId", "entrezId"),
                y = colnames(mcols(object))
            )
        )
        df <- mcols(object)
        colnames(df)[colnames(df) == "geneId"] <- "ensemblId"
        metadata(df) <- metadata(object)
        .makeEnsembl2Entrez(
            object = df,
            format = match.arg(format)
        )
    }

formals(`Ensembl2Entrez,GenomicRanges`)[["format"]] <- # nolint
    formals(.makeEnsembl2Entrez)[["format"]]



#' @rdname Ensembl2Entrez
#' @export
setMethod(
    f = "Ensembl2Entrez",
    signature = signature(object = "GenomicRanges"),
    definition = `Ensembl2Entrez,GenomicRanges`
)

#' @rdname Ensembl2Entrez
#' @export
setMethod(
    f = "Ensembl2Entrez",
    signature = signature(object = "character"),
    definition = `Ensembl2Entrez,character`
)
