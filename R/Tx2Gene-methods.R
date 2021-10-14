#' @name Tx2Gene
#' @inherit Tx2Gene-class title description return
#' @note Updated 2021-08-10.
#'
#' @inheritParams AcidRoxygen::params
#' @inheritParams params
#' @param ... Arguments pass through to `DataFrame` method.
#'
#' @seealso [makeTx2Gene()].
#'
#' @examples
#' ## DataFrame ====
#' df <- DataFrame(
#'     "txId" = c(
#'         "ENST00000635602.1",
#'         "ENST00000635506.1"
#'     ),
#'     "geneId" = c(
#'         "ENSG00000283061.1",
#'         "ENSG00000283061.1"
#'     )
#' )
#' t2g <- Tx2Gene(df)
#' print(t2g)
NULL



## Updated 2021-10-13.
`Tx2Gene,DFrame` <-  # nolint
    function(
        object,
        quiet = FALSE
    ) {
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
            hasNoDuplicates(object[["txId"]]),
            msg = "Failed to generate Tx2Gene object."
        )
        object <- object[order(object), , drop = FALSE]
        metadata(object) <- meta
        new(Class = "Tx2Gene", object)
    }



## Updated 2021-08-09.
`Tx2Gene,GRanges` <-  # nolint
    function(object, ...) {
        df <- as(object, "DFrame")
        metadata(df) <- metadata(object)
        Tx2Gene(df, ...)
    }



## Updated 2021-08-09.
`Tx2Gene,GRangesList` <-  # nolint
    function(object, ...) {
        x <- as(object, "GRangesList")
        x <- unname(x)
        gr <- unlist(x)
        assert(is(gr, "GRanges"))
        Tx2Gene(gr, ...)
    }



## Updated 2021-08-09.
`Tx2Gene,data.frame` <-  # nolint
    function(object, ...) {
        assert(identical(ncol(object), 2L))
        colnames(object) <- c("txId", "geneId")
        object <- as(object, "DFrame")
        Tx2Gene(object, ...)
    }



## Updated 2021-08-09.
`Tx2Gene,matrix` <-  # nolint
    function(object, ...) {
        assert(
            is.character(object),
            identical(ncol(object), 2L)
        )
        object <- as.data.frame(object, stringsAsFactors = FALSE)
        Tx2Gene(object, ...)
    }



#' @rdname Tx2Gene
#' @export
setMethod(
    f = "Tx2Gene",
    signature = signature(object = "DFrame"),
    definition = `Tx2Gene,DFrame`
)

#' @rdname Tx2Gene
#' @export
setMethod(
    f = "Tx2Gene",
    signature = signature(object = "GRanges"),
    definition = `Tx2Gene,GRanges`
)

#' @rdname Tx2Gene
#' @export
setMethod(
    f = "Tx2Gene",
    signature = signature(object = "GRangesList"),
    definition = `Tx2Gene,GRangesList`
)

#' @rdname Tx2Gene
#' @export
setMethod(
    f = "Tx2Gene",
    signature = signature(object = "data.frame"),
    definition = `Tx2Gene,data.frame`
)

#' @rdname Tx2Gene
#' @export
setMethod(
    f = "Tx2Gene",
    signature = signature(object = "matrix"),
    definition = `Tx2Gene,matrix`
)
