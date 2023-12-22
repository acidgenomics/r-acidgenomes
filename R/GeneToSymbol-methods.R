#' @name GeneToSymbol
#' @inherit AcidGenerics::GeneToSymbol description return title
#' @note Updated 2023-12-21.
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
#' @seealso [makeGeneToSymbol()].
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
#' x <- GeneToSymbol(df)
#' print(x)
#'
#' ## GRanges ====
#' object <- GRanges
#' x <- GeneToSymbol(object)
NULL



## Updated 2023-12-21.
`GeneToSymbol,DFrame` <- # nolint
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
        object <- as(object, "DFrame")
        cols <- c("geneId", "geneName")
        if (!isSubset(cols, colnames(object))) {
            colnames(object) <- camelCase(colnames(object), strict = TRUE)
        }
        assert(
            isSubset(cols, colnames(object)),
            hasRows(object),
            !anyNA(object[["geneId"]])
        )
        object <- object[, cols, drop = FALSE]
        object <- decode(object)
        object <- unique(object)
        keep <- complete.cases(object)
        if (!all(keep)) {
            ## e.g. applies to Ensembl Mus musculus GRCm39 104.
            meta[["dropped"]] <- object[["geneId"]][!keep]
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
        assert(
            hasRows(object),
            hasNoDuplicates(object[["geneId"]])
        )
        i <- order(object)
        object <- object[i, , drop = FALSE]
        dupes <- sort(dupes(object[["geneName"]]))
        if (hasLength(dupes)) {
            meta[["dupes"]] <- dupes
        }
        switch(
            EXPR = format,
            "makeUnique" = {
                object[["geneName"]] <- make.unique(object[["geneName"]])
            },
            "1:1" = {
                i <- order(object[["geneName"]], object[["geneId"]])
                object <- object[i, , drop = FALSE]
                i <- !duplicated(object[["geneName"]])
                object <- object[i, , drop = FALSE]
            },
            "unmodified" = {
                if (isFALSE(quiet) && hasLength(dupes)) {
                    alertInfo(sprintf(
                        "Returning unmodified with %d duplicate gene %s: %s.",
                        length(dupes),
                        ngettext(
                            n = length(dupes),
                            msg1 = "symbol",
                            msg2 = "symbols"
                        ),
                        toInlineString(dupes)
                    ))
                }
            }
        )
        i <- order(object)
        object <- object[i, , drop = FALSE]
        meta[["date"]] <- Sys.Date()
        meta[["packageVersion"]] <- .pkgVersion
        meta <- meta[sort(names(meta))]
        metadata(object) <- meta
        new(Class = "GeneToSymbol", object)
    }



## Updated 2023-04-26.
`GeneToSymbol,GRanges` <- # nolint
    function(object, ...) {
        df <- as(object, "DFrame")
        metadata(df) <- metadata(object)
        GeneToSymbol(df, ...)
    }



## Updated 2023-12-21.
`GeneToSymbol,GRangesList` <- # nolint
    function(object, ...) {
        gr <- unlist(x = object, recursive = FALSE, use.names = TRUE)
        GeneToSymbol(gr, ...)
    }



#' @rdname GeneToSymbol
#' @export
setMethod(
    f = "GeneToSymbol",
    signature = signature(object = "DFrame"),
    definition = `GeneToSymbol,DFrame`
)

#' @rdname GeneToSymbol
#' @export
setMethod(
    f = "GeneToSymbol",
    signature = signature(object = "GRanges"),
    definition = `GeneToSymbol,GRanges`
)

#' @rdname GeneToSymbol
#' @export
setMethod(
    f = "GeneToSymbol",
    signature = signature(object = "GRangesList"),
    definition = `GeneToSymbol,GRangesList`
)
