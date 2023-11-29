#' @name TxToGene
#' @inherit AcidGenerics::TxToGene description return title
#' @note Updated 2023-11-29.
#'
#' @inheritParams AcidRoxygen::params
#' @inheritParams params
#' @param ... Arguments pass through to `DFrame` method.
#'
#' @seealso [makeTxToGene()].
#'
#' @examples
#' ## DFrame ====
#' df <- S4Vectors::DataFrame(
#'     "txId" = c(
#'         "ENST00000635602.1",
#'         "ENST00000635506.1"
#'     ),
#'     "geneId" = c(
#'         "ENSG00000283061.1",
#'         "ENSG00000283061.1"
#'     )
#' )
#' object <- TxToGene(df)
#' print(object)
NULL



## Updated 2023-11-29.
`TxToGene,DFrame` <- # nolint
    function(object,
             quiet = FALSE) {
        assert(
            hasColnames(object),
            hasRows(object),
            isFlag(quiet)
        )
        meta <- metadata(object)
        cols <- c("txId", "geneId")
        ## Rename legacy columns, if necessary.
        if (!isSubset(cols, colnames(object))) {
            colnames(object) <- camelCase(colnames(object), strict = TRUE)
            colnames(object) <- gsub(
                pattern = "^transcript",
                replacement = "tx",
                x = colnames(object)
            )
        }
        assert(
            isSubset(cols, colnames(object)),
            hasRows(object)
        )
        object <- object[, cols, drop = FALSE]
        object <- decode(object)
        assert(allAreAtomic(object))
        object <- unique(object)
        ## Allowing sanitization of messy input with incomplete elements.
        keep <- complete.cases(object)
        if (!all(keep)) {
            meta[["dropped"]] <- which(!keep)
            if (isFALSE(quiet)) {
                n <- sum(!keep)
                alertWarning(sprintf(
                    "Dropping %d %s without transcript-to-gene mapping.",
                    n,
                    ngettext(
                        n = n,
                        msg1 = "element",
                        msg2 = "elements"
                    )
                ))
            }
            object <- object[keep, , drop = FALSE]
        }
        assert(
            hasRows(object),
            all(complete.cases(object)),
            hasNoDuplicates(object[["txId"]])
        )
        object <- object[order(object), , drop = FALSE]
        meta[["date"]] <- Sys.Date()
        meta[["packageVersion"]] <- .pkgVersion
        meta <- meta[sort(names(meta))]
        metadata(object) <- meta
        new(Class = "TxToGene", object)
    }



## Updated 2021-08-09.
`TxToGene,GRanges` <- # nolint
    function(object, ...) {
        df <- as(object, "DFrame")
        metadata(df) <- metadata(object)
        TxToGene(df, ...)
    }



## Updated 2023-04-26.
`TxToGene,GRangesList` <- # nolint
    function(object, ...) {
        grl <- as(object, "GRangesList")
        gr <- unlist(grl, recursive = FALSE, use.names = TRUE)
        assert(is(gr, "GRanges"))
        TxToGene(gr, ...)
    }



## Updated 2021-08-09.
`TxToGene,data.frame` <- # nolint
    function(object, ...) {
        assert(identical(ncol(object), 2L))
        colnames(object) <- c("txId", "geneId")
        object <- as(object, "DFrame")
        TxToGene(object, ...)
    }



## Updated 2021-08-09.
`TxToGene,matrix` <- # nolint
    function(object, ...) {
        assert(
            is.character(object),
            identical(ncol(object), 2L)
        )
        object <- as.data.frame(object, stringsAsFactors = FALSE)
        TxToGene(object, ...)
    }



#' @rdname TxToGene
#' @export
setMethod(
    f = "TxToGene",
    signature = signature(object = "DFrame"),
    definition = `TxToGene,DFrame`
)

#' @rdname TxToGene
#' @export
setMethod(
    f = "TxToGene",
    signature = signature(object = "GRanges"),
    definition = `TxToGene,GRanges`
)

#' @rdname TxToGene
#' @export
setMethod(
    f = "TxToGene",
    signature = signature(object = "GRangesList"),
    definition = `TxToGene,GRangesList`
)

#' @rdname TxToGene
#' @export
setMethod(
    f = "TxToGene",
    signature = signature(object = "data.frame"),
    definition = `TxToGene,data.frame`
)

#' @rdname TxToGene
#' @export
setMethod(
    f = "TxToGene",
    signature = signature(object = "matrix"),
    definition = `TxToGene,matrix`
)
