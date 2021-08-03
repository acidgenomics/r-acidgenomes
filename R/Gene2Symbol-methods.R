## FIXME Need to ensure we handle NA values correctly here.



#' @name Gene2Symbol
#' @inherit Gene2Symbol-class title description return
#'
#' @note For some organisms, gene names and gene symbols do not map 1:1 (e.g.
#'   *Homo sapiens* and *Mus musculus*). Refer to the `format` argument here in
#'   the documentation for approaches that deal with this issue.
#' @note For the `format` argument, note that "long" was used instead of
#'   "unmodified" prior to v0.10.10.
#' @note Updated 2021-06-09.
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
#' data(GRanges, package = "AcidTest")
#'
#' ## DataFrame ====
#' df <- DataFrame(
#'     "geneId" = c(
#'         "ENSG00000228572.7",
#'         "ENSG00000182378.14"
#'     ),
#'     "geneName" = c(
#'         "AL954722.1",
#'         "PLCXD1"
#'     )
#' )
#' x <- Gene2Symbol(df)
#' print(x)
#'
#' ## GRanges ====
#' object <- GRanges
#' x <- Gene2Symbol(object)
NULL



## Updated 2021-06-09.
`Gene2Symbol,DataFrame` <-  # nolint
    function(object, format = c("makeUnique", "unmodified", "1:1")) {
        format <- match.arg(format)
        cols <- c("geneId", "geneName")
        if (!isSubset(cols, colnames(object))) {
            colnames(object) <- camelCase(colnames(object), strict = TRUE)
        }
        assert(
            hasRows(object),
            isSubset(cols, colnames(object))
        )
        df <- decode(object[, cols])
        ## Allow coercion of integer gene identifiers (e.g. NCBI Entrez).
        if (is.integer(df[[cols[[1L]]]])) {
            df[[cols[[1L]]]] <- as.character(df[[cols[[1L]]]])
        }
        ## Inform the user about how many symbols multi-map.
        ## Note that `duplicated()` doesn't work on Rle, so we have to decode.
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
        switch(
            EXPR = format,
            "makeUnique" = {
                ## Returning 1:1 mappings with renamed gene symbols. This is the
                ## default, and including a message is too noisy, since it is
                ## used heavily in other functions.
                df[["geneName"]] <- make.unique(df[["geneName"]])
            },
            "unmodified" = {
                alertWarning(paste(
                    "Returning with unmodified gene symbols",
                    "{.emph (may contain duplicates)}."
                ))
            },
            "1:1" = {
                alert(paste(
                    "Returning 1:1 mappings using oldest",
                    "gene identifier per symbol."
                ))
                x <- split(x = df, f = df[[1L]])
                assert(is(x, "SplitDataFrameList"))
                x <- lapply(
                    X = x[, 2L],
                    FUN = function(x) {
                        sort(x = x, decreasing = FALSE, na.last = TRUE)[[1L]]
                    }
                )
                x <- unlist(x, recursive = FALSE, use.names = TRUE)
                df <- DataFrame("a" = names(x), "b" = x)
                rownames(df) <- df[[1L]]
                colnames(df) <- cols
            }
        )
        assert(is(df, "DataFrame"))
        metadata(df) <- metadata(object)
        metadata(df)[["format"]] <- format
        metadata(df)[c("call", "synonyms")] <- NULL
        new(Class = "Gene2Symbol", df)
    }



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
    signature = signature("DataFrame"),
    definition = `Gene2Symbol,DataFrame`
)

#' @rdname Gene2Symbol
#' @export
setMethod(
    f = "Gene2Symbol",
    signature = signature("GRanges"),
    definition = `Gene2Symbol,GRanges`
)
