#' @name Gene2Symbol
#' @inherit AcidGenerics::Gene2Symbol description return title
#' @note Updated 2023-09-16.
#'
#' @details
#' For some organisms, gene identifiers and gene names do not map 1:1 (e.g.
#' *Homo sapiens* and *Mus musculus*). Refer to the `format` argument
#' here in the documentation for approaches that deal with this issue.
#'
#' @inheritParams AcidRoxygen::params
#' @inheritParams params
#'
#' @param format `character(1)`.
#' Formatting method to apply:
#'
#' - `"makeUnique"`: *Recommended.* Apply `make.unique` to the `geneName`
#' column. Gene names are made unique, while the identifiers remain
#' unmodified. `NA` gene names will be renamed to `"unannotated"`.
#' - `"1:1"`: For gene names that map to multiple gene identifiers, select
#' only the first annotated gene identifier. Incomplete elements with
#' `NA` gene name will be removed will be removed with an internal
#' `complete.cases` call.
#' - `"unmodified"`: Return `geneId` and `geneName` columns unmodified, in
#' long format. Incomplete elements with `NA` gene name will be removed
#' with an internal `complete.cases` call.
#'
#' @param ... Arguments pass through to `DFrame` method.
#'
#' @seealso [makeGene2Symbol()].
#'
#' @examples
#' data(GRanges, package = "AcidTest")
#'
#' ## DFrame ====
#' df <- S4Vectors::DataFrame(
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



## Updated 2023-09-16.
`Gene2Symbol,DFrame` <- # nolint
    function(object,
             format = c("makeUnique", "1:1", "unmodified"),
             quiet = FALSE) {
        assert(
            hasColnames(object),
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
        object <- as(object, "DFrame")
        object <- object[, cols, drop = FALSE]
        object <- decode(object)
        assert(allAreAtomic(object))
        object <- unique(object)
        ## Remove incomplete elements, except for "makeUnique" mode.
        if (!identical(format, "makeUnique")) {
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
        assert(hasRows(object))
        ## Enforce coercion of integer gene identifiers (e.g. NCBI).
        if (is.integer(object[["geneId"]])) {
            object[["geneId"]] <- as.character(object[["geneId"]])
        }
        switch(
            EXPR = format,
            "makeUnique" = {
                ## Ensure we arrange by gene identifier prior to sanization.
                object <- object[order(object), , drop = FALSE]
                ## Replace any "NA" gene names with "unannotated".
                if (anyNA(object[["geneName"]])) {
                    object[["geneName"]][
                        which(is.na(object[["geneName"]]))
                    ] <- "unannotated"
                }
                ## Inform the user about how many symbols multi-map.
                ## Don't count "unannotated" genes as true duplicates here.
                ## Note that `duplicated()` doesn't work on Rle.
                dupes <- duplicated(object[["geneName"]])
                if (any(dupes)) {
                    dupes <- setdiff(
                        x = sort(unique(object[["geneName"]][dupes])),
                        y = "unannotated"
                    )
                    if (hasLength(dupes)) {
                        meta[["dupes"]] <- dupes
                        if (isFALSE(quiet)) {
                            n <- length(dupes)
                            alertInfo(sprintf(
                                "%d non-unique gene %s detected: %s.",
                                n,
                                ngettext(
                                    n = n,
                                    msg1 = "symbol",
                                    msg2 = "symbols"
                                ),
                                toInlineString(dupes, n = 10L, class = "val")
                            ))
                        }
                    }
                }
                object[["geneName"]] <- make.unique(object[["geneName"]])
            },
            "1:1" = {
                assert(all(complete.cases(object)))
                alert(paste(
                    "Returning 1:1 mappings using oldest",
                    "gene identifier per symbol."
                ))
                x <- split(x = object, f = object[["geneName"]])
                assert(is(x, "SplitDFrameList"))
                x <- SplitDataFrameList(lapply(
                    X = x,
                    FUN = function(x) {
                        idx <- order(
                            x = x[["geneId"]],
                            decreasing = FALSE,
                            na.last = TRUE
                        )
                        head(x[idx, , drop = FALSE], n = 1L)
                    }
                ))
                object <- do.call(what = rbind, args = x)
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
            is(object, "DFrame"),
            all(complete.cases(object)),
            hasNoDuplicates(object[["geneId"]]),
            msg = "Failed to generate Gene2Symbol object."
        )
        object <- object[order(object), , drop = FALSE]
        meta[["date"]] <- Sys.Date()
        meta[["packageVersion"]] <- .pkgVersion
        meta <- meta[sort(names(meta))]
        metadata(object) <- meta
        new(Class = "Gene2Symbol", object)
    }



## Updated 2023-04-26.
`Gene2Symbol,GRanges` <- # nolint
    function(object, ...) {
        df <- as(object, "DFrame")
        metadata(df) <- metadata(object)
        Gene2Symbol(df, ...)
    }



#' @rdname Gene2Symbol
#' @export
setMethod(
    f = "Gene2Symbol",
    signature = signature(object = "DFrame"),
    definition = `Gene2Symbol,DFrame`
)

#' @rdname Gene2Symbol
#' @export
setMethod(
    f = "Gene2Symbol",
    signature = signature(object = "GRanges"),
    definition = `Gene2Symbol,GRanges`
)
