#' @inherit Ensembl2Entrez-class title description return
#' @name Ensembl2Entrez
#' @note Updated 2020-10-01.
#'
#' @inheritParams AcidRoxygen::params
#' @param format `character(1)`.
#'   Formatting method to apply:
#'
#'   - `"1:1"`: *Recommended.*
#'       Return with 1:1 mappings. For Ensembl genes that don't map 1:1 with
#'       Entrez, pick the oldest Entrez identifier. Genes that don't map to
#'       Entrez will contain `NA` in `entrezId` column.
#'   - `"long"`:
#'       Return 1:many in long format.
#'
#' @examples
#' data(RangedSummarizedExperiment, package = "AcidTest")
#' rse <- RangedSummarizedExperiment
#'
#' ## character ====
#' ## Ensembl-to-Entrez.
#' genes <- c("ENSG00000000003", "ENSG00000000005")
#' x <- Ensembl2Entrez(genes)
#' print(x)
#'
#' ## integer ====
#' ## Entrez-to-Ensembl.
#' genes <- c(1L, 2L)
#' x <- Entrez2Ensembl(genes)
#' print(x)
#'
#' ## SummarizedExperiment ====
#' x <- Ensembl2Entrez(rse)
#' print(x)
NULL



#' Make an Ensembl2Entrez (or Entrez2Ensembl) object
#'
#' @note Updated 2019-08-16.
#' @noRd
.makeEnsembl2Entrez <-
    function(
        object,
        format = c("1:1", "long"),
        return = c("Ensembl2Entrez", "Entrez2Ensembl")
    ) {
        format <- match.arg(format)
        return <- match.arg(return)
        cols <- switch(
            EXPR = return,
            "Ensembl2Entrez" = c("ensembl", "entrez"),
            "Entrez2Ensembl" = c("entrez", "ensembl")
        )
        assert(
            is(object, "DataFrame"),
            hasRows(object),
            isSubset(cols, colnames(object))
        )
        df <- object[, cols]
        rownames(df) <- NULL
        df <- decode(expand(df))
        assert(
            is.character(df[["ensembl"]]),
            is.integer(df[["entrez"]])
        )
        if (identical(format, "1:1")) {
            split <- split(x = df, f = df[[1L]])
            unique <- unlist(
                x = bplapply(
                    X = split[, 2L],
                    FUN = function(x) {
                        if (all(is.na(x))) {
                            NA
                        } else {
                            head(sort(x), n = 1L)
                        }
                    }
                ),
                recursive = FALSE,
                use.names = FALSE
            )
            df <- DataFrame("a" = names(split), "b" = unique)
            ## This will check the return values for duplicates.
            ## > assert(all(bapply(X = df, FUN = hasNoDuplicates)))
            rownames(df) <- df[[1L]]
            colnames(df) <- cols
        }
        df[["entrez"]] <- as.integer(df[["entrez"]])
        df <- df[complete.cases(df), ]
        assert(hasRows(df))
        metadata(df) <- metadata(object)
        metadata(df)[["format"]] <- format
        new(Class = return, df)
    }



## Updated 2020-10-01.
`Ensembl2Entrez,character` <-  # nolint
    function(
        object,
        organism = NULL,
        format = c("1:1", "long")
    ) {
        if (is.null(organism)) {
            organism <- detectOrganism(object)
        }
        .makeEnsembl2Entrez(
            object = .getEnsembl2EntrezFromOrgDb(
                keys = object,
                keytype = "ENSEMBL",
                columns = "ENTREZID",
                organism = organism
            ),
            format = match.arg(format)
        )
    }

#' @rdname Ensembl2Entrez
#' @export
setMethod(
    f = "Ensembl2Entrez",
    signature = signature("character"),
    definition = `Ensembl2Entrez,character`
)



## Updated 2020-10-01.
`Entrez2Ensembl,integer` <-  # nolint
    function(
        object,
        organism = "Homo sapiens",
        format = c("1:1", "long")
    ) {
        .makeEnsembl2Entrez(
            object = .getEnsembl2EntrezFromOrgDb(
                keys = as.character(object),
                keytype = "ENTREZID",
                columns = "ENSEMBL",
                organism = organism
            ),
            format = match.arg(format),
            return = "Entrez2Ensembl"
        )
    }

#' @rdname Ensembl2Entrez
#' @export
setMethod(
    f = "Entrez2Ensembl",
    signature = signature("integer"),
    definition = `Entrez2Ensembl,integer`
)



## Updated 2020-10-01.
`Ensembl2Entrez,GRanges` <-  # nolint
    function(
        object,
        format = c("1:1", "long")
    ) {
        df <- as(object, "DataFrame")
        colnames(df)[colnames(df) == "geneId"] <- "ensembl"
        colnames(df)[colnames(df) == "entrezId"] <- "entrez"
        metadata(df) <- metadata(object)
        .makeEnsembl2Entrez(
            object = df,
            format = match.arg(format)
        )
    }



#' @rdname Ensembl2Entrez
#' @export
setMethod(
    f = "Ensembl2Entrez",
    signature = signature("GRanges"),
    definition = `Ensembl2Entrez,GRanges`
)



## Updated 2020-10-01.
`Ensembl2Entrez,RangedSummarizedExperiment` <-  # nolint
    function(object, ...) {
        Ensembl2Entrez(rowRanges(object), ...)
    }



#' @rdname Ensembl2Entrez
#' @export
setMethod(
    f = "Ensembl2Entrez",
    signature = signature("RangedSummarizedExperiment"),
    definition = `Ensembl2Entrez,RangedSummarizedExperiment`
)
