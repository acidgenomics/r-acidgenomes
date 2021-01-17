#' @inherit Gene2Symbol-class title description return
#' @name Gene2Symbol
#'
#' @note For some organisms, gene names and gene symbols do not map 1:1 (e.g.
#'   *Homo sapiens* and *Mus musculus*). Refer to the `format` argument here in
#'   the documentation for approaches that deal with this issue.
#' @note For the `format` argument, note that "long" was used instead of
#'   "unmodified" prior to v0.10.10.
#' @note Updated 2021-01-17.
#'
#' @inheritParams AcidRoxygen::params
#' @param format `character(1)`.
#'   Formatting method to apply:
#'
#'   - `"makeUnique"`: *Recommended.* Apply [`make.unique()`][base::make.unique]
#'     to the `geneName` column. Gene symbols are made unique, while the gene
#'     identifiers remain unmodified.
#'   - `"unmodified"`: Return `geneId` and `geneName` columns unmodified, in
#'     long format.
#'   - `"1:1"`: For gene symbols that map to multiple gene identifiers, select
#"     only the first annotated gene identifier.
#'
#' @seealso [makeGene2Symbol()].
#'
#' @examples
#' data(RangedSummarizedExperiment, package = "AcidTest")
#' rse <- RangedSummarizedExperiment
#'
#' ## SummarizedExperiment ====
#' x <- Gene2Symbol(rse)
#' print(x)
NULL



## Updated 2021-01-17.
`Gene2Symbol,DataFrame` <-  # nolint
    function(object, format = c("makeUnique", "unmodified", "1:1")) {
        assert(hasRows(object))
        format <- match.arg(format)
        colnames(object) <- camelCase(colnames(object), strict = TRUE)
        cols <- c("geneId", "geneName")
        if (!all(cols %in% colnames(object))) {
            stop(sprintf(
                "Object does not contain gene-to-symbol mappings: %s.",
                toString(cols)
            ))
        }
        df <- DataFrame(
            "geneId" = as.character(decode(object[["geneId"]])),
            "geneName" = as.character(decode(object[["geneName"]])),
            row.names = rownames(object)
        )
        ## Inform the user about how many symbols multi-map.
        ## Note that `duplicated` doesn't work on Rle, so we have to coerce
        ## columns to character first (see `as_tibble` call above).
        duplicated <- duplicated(df[["geneName"]])
        if (any(duplicated)) {
            dupes <- unique(df[["geneName"]][duplicated])
            alertInfo(sprintf(
                "%d non-unique gene %s detected.",
                length(dupes),
                ngettext(
                    n = length(dupes),
                    msg1 = "symbol",
                    msg2 = "symbols"
                )
            ))
        }
        ## Return mode.
        if (format == "makeUnique") {
            ## Returning 1:1 mappings with renamed gene symbols.
            ## This is the default, and including a message is too noisy, since
            ## it is used heavily in other functions.
            df[["geneName"]] <- make.unique(df[["geneName"]])
        } else if (format == "unmodified") {
            alertWarning(paste(
                "Returning with unmodified gene symbols",
                "{.emph (may contain duplicates)}."
            ))
        } else if (format == "1:1") {
            alert(paste(
                "Returning 1:1 mappings using oldest",
                "gene identifier per symbol."
            ))
            x <- split(df, f = df[["geneName"]])
            x <- bplapply(
                X = x,
                FUN = function(x) {
                    x <- x[order(x[["geneId"]]), , drop = FALSE]
                    x <- head(x, n = 1L)
                    x
                }
            )
            x <- DataFrameList(x)
            x <- unlist(x, recursive = FALSE, use.names = FALSE)
            df <- x
            assert(is(df, "DataFrame"))
        }
        metadata(df) <- .slotGenomeMetadata(object)
        metadata(df)[["format"]] <- format
        new(Class = "Gene2Symbol", df)
    }



#' @rdname Gene2Symbol
#' @export
setMethod(
    f = "Gene2Symbol",
    signature = signature("DataFrame"),
    definition = `Gene2Symbol,DataFrame`
)



## Updated 2019-07-22.
`Gene2Symbol,GRanges` <-  # nolint
    function(object, format) {
        df <- as(object, "DataFrame")
        df <- unique(df)
        metadata(df) <- metadata(object)
        do.call(what = Gene2Symbol, args = list(object = df, format = format))
    }

formals(`Gene2Symbol,GRanges`) <- formals(`Gene2Symbol,DataFrame`)



#' @rdname Gene2Symbol
#' @export
setMethod(
    f = "Gene2Symbol",
    signature = signature("GRanges"),
    definition = `Gene2Symbol,GRanges`
)



## Updated 2019-07-22.
`Gene2Symbol,SummarizedExperiment` <-  # nolint
    function(object, format) {
        object <- as.SummarizedExperiment(object)
        df <- rowData(object)
        rownames(df) <- rownames(object)
        do.call(what = Gene2Symbol, args = list(object = df, format = format))
    }

formals(`Gene2Symbol,SummarizedExperiment`) <- formals(`Gene2Symbol,DataFrame`)



#' @rdname Gene2Symbol
#' @export
setMethod(
    f = "Gene2Symbol",
    signature = signature("SummarizedExperiment"),
    definition = `Gene2Symbol,SummarizedExperiment`
)
