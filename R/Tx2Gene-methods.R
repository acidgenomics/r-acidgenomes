#' @name Tx2Gene
#' @inherit Tx2Gene-class title description return
#'
#' @note No attempt is made to arrange the rows by transcript identifier.
#' @note Updated 2021-08-09.
#'
#' @inheritParams AcidRoxygen::params
#' @param removeNA `logical(1)`.
#'   Remove elements that don't contain clear mappings (e.g. gene identifier
#'   is `NA`). Calls `complete.cases` internally.
#'   If set `FALSE`, the function will intentionally error if any incomplete
#'   cases are detected.
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



## Updated 2021-08-09.
`Tx2Gene,data.frame` <-  # nolint
    function(object, ...) {
        assert(identical(ncol(object), 2L))
        colnames(object) <- c("txId", "geneId")
        object <- as(object, "DataFrame")
        Tx2Gene(object, ...)
    }



## Updated 2021-08-09.
`Tx2Gene,DataFrame` <-  # nolint
    function(
        object,
        removeNA = FALSE
    ) {
        assert(hasColnames(object))
        meta <- metadata(object)
        cols <- c("txId", "geneId")
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
        if (isTRUE(removeNA)) {
            keep <- complete.cases(object)
            if (!all(keep)) {
                meta[["dropped"]] <- which(!keep)
                n <- sum(!keep)
                alertWarning(sprintf(
                    "Dropping %d %s without transcript-to-gene mapping.",
                    n,
                    ngettext(
                        n = n,
                        msg1 = "identifier",
                        msg2 = "identifiers"
                    )
                ))
                object <- object[keep, , drop = FALSE]
            }
        }
        assert(
            hasRows(object),
            hasNoDuplicates(object[[cols[[1L]]]]),
            all(complete.cases(object))
        )
        object <- object[order(object), , drop = FALSE]
        metadata(object) <- meta
        new(Class = "Tx2Gene", object)
    }



## Updated 2021-08-09.
`Tx2Gene,GRanges` <-  # nolint
    function(object, ...) {
        df <- as(object, "DataFrame")
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



#' @rdname Tx2Gene
#' @export
setMethod(
    f = "Tx2Gene",
    signature = signature("DataFrame"),
    definition = `Tx2Gene,DataFrame`
)

#' @rdname Tx2Gene
#' @export
setMethod(
    f = "Tx2Gene",
    signature = signature("GRanges"),
    definition = `Tx2Gene,GRanges`
)

#' @rdname Tx2Gene
#' @export
setMethod(
    f = "Tx2Gene",
    signature = signature("GRangesList"),
    definition = `Tx2Gene,GRangesList`
)

#' @rdname Tx2Gene
#' @export
setMethod(
    f = "Tx2Gene",
    signature = signature("data.frame"),
    definition = `Tx2Gene,data.frame`
)

#' @rdname Tx2Gene
#' @export
setMethod(
    f = "Tx2Gene",
    signature = signature("matrix"),
    definition = `Tx2Gene,matrix`
)
