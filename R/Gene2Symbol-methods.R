## FIXME 1:1 isn't returning unique for gene symbols.
## FIXME Need to think about this better for rRNA genes, etc.
##
## Use "rse.rds" example on desktop to better check this code. Need to
## document more clearly, and ensure we recommend makeUnique by default.
## FIXME makeUnique approach should convert NA_character_ to "NA", correct?
## FIXME Dealing with this is such a pain...

## FIXME Need to rework the "unmodified" handling here.



#' @name Gene2Symbol
#' @inherit Gene2Symbol-class title description return
#'
#' @note For some organisms, gene names and gene symbols do not map 1:1 (e.g.
#'   *Homo sapiens* and *Mus musculus*). Refer to the `format` argument here in
#'   the documentation for approaches that deal with this issue.
#' @note For the `format` argument, note that "long" was used instead of
#'   "unmodified" prior to v0.10.10.
#' @note Updated 2021-08-09.
#'
#' @inheritParams AcidRoxygen::params
#' @inheritParams params
#' @param format `character(1)`.
#'   Formatting method to apply:
#'
#'   - `"makeUnique"`: *Recommended.* Apply `make.unique` to the `geneName`
#'     column. Gene symbols are made unique, while the identifiers remain
#'     unmodified.
#'   - `"unmodified"`: Return `geneId` and `geneName` columns unmodified, in
#'     long format.
#'   - `"1:1"`: For gene symbols that map to multiple gene identifiers, select
#'     only the first annotated gene identifier.
#' @param ... Arguments pass through to `DataFrame` method.
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



## NOTE "unmodified" argument is used by `AcidExperiment::mapGenes`.
## Updated 2021-08-10.
`Gene2Symbol,DataFrame` <-  # nolint
    function(
        object,
        format = c("makeUnique", "1:1", "unmodified"),
        completeCases = TRUE,
        quiet = FALSE
    ) {
        assert(
            hasColnames(object),
            isFlag(completeCases),
            isFlag(quiet)
        )
        format <- match.arg(format)
        meta <- metadata(object)
        meta[["format"]] <- format
        cols <- c("geneId", "geneName")
        if (!isSubset(cols, colnames(object))) {
            colnames(object) <- camelCase(colnames(object), strict = TRUE)
        }
        assert(
            isSubset(cols, colnames(object)),
            hasRows(object)
        )
        object <- as(object, "DataFrame")
        object <- object[, cols, drop = FALSE]
        object <- decode(object)
        assert(allAreAtomic(object))
        object <- unique(object)
        if (isTRUE(completeCases)) {
            keep <- complete.cases(object)
            if (!all(keep)) {
                ## e.g. applies to Ensembl Mus musculus GRCm39 104.
                meta[["dropped"]] <- which(!keep)
                if (isFALSE(quiet)) {
                    n <- sum(!keep)
                    alertWarning(sprintf(
                        "Dropping %d %s without defined gene symbol.",
                        n,
                        ngettext(
                            n = n,
                            msg1 = "identifier",
                            msg2 = "identifiers"
                        )
                    ))
                }
                object <- object[keep, , drop = FALSE]
            }
        }

        ## FIXME Don't use numeric columns here below...
        ## FIXME Also change this convention in Tx2Gene code...

        assert(hasRows(object))
        ## Enforce coercion of integer gene identifiers (e.g. NCBI Entrez).
        if (is.integer(object[[cols[[1L]]]])) {
            object[[cols[[1L]]]] <- as.character(object[[cols[[1L]]]])
        }
        switch(
            EXPR = format,
            "makeUnique" = {
                ## Inform the user about how many symbols multi-map.
                ## Note that `duplicated()` doesn't work on Rle.
                dupes <- duplicated(object[["geneName"]])
                if (any(dupes)) {
                    dupes <- sort(unique(object[["geneName"]][dupes]))
                    meta[["dupes"]] <- dupes
                    if (isFALSE(quiet)) {
                        n <- length(dupes)
                        alertInfo(sprintf(
                            "%d non-unique gene %s detected.",
                            n,
                            ngettext(
                                n = n,
                                msg1 = "symbol",
                                msg2 = "symbols"
                            )
                        ))
                    }
                }
                object[["geneName"]] <- make.unique(object[["geneName"]])
            },
            "1:1" = {
                ## FIXME Handle the dupes quietly here...
                ## FIXME Don't think this is working quite right for 1:1
                ## mappings.
                ## FIXME This isn't working the way we want...
                ## FIXME Check handling of ribosomal genes.
                ## FIXME Need to split by the gene symbol
                alert(paste(
                    "Returning 1:1 mappings using oldest",
                    "gene identifier per symbol."
                ))
                x <- split(x = object, f = object[[1L]])
                assert(is(x, "SplitDataFrameList"))
                x <- lapply(
                    X = x[, 2L],
                    FUN = function(x) {
                        sort(x = x, decreasing = FALSE, na.last = TRUE)[[1L]]
                    }
                )
                x <- unlist(x, recursive = FALSE, use.names = TRUE)
                object <- DataFrame(names(x), x)
                dimnames(object) <- list(object[[1L]], cols)
            },
            "unmodified" = {
                if (isFALSE(quiet)) {
                    alertInfo(paste(
                        "Returning with unmodified gene symbols",
                        "{.emph (may contain duplicates)}."
                    ))
                }
            }
        )
        assert(
            is(object, "DataFrame"),
            all(complete.cases(object)),
            hasNoDuplicates(object[[cols[[1L]]]]),
            msg = "Failed to generate Gene2Symbol object."
        )
        object <- object[order(object), , drop = FALSE]
        metadata(object) <- meta
        new(Class = "Gene2Symbol", object)
    }



## Updated 2021-08-09.
`Gene2Symbol,GRanges` <-  # nolint
    function(object, ...) {
        df <- as(object, "DataFrame")
        metadata(df) <- metadata(object)
        Gene2Symbol(df, ...)
    }



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
