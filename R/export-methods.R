#' @name export
#' @inherit pipette::export description params return title
#' @note Updated 2022-09-13.
#'
#' @details
#' The `Tx2Gene` method automatically disables writing of column names, which
#' is the intended input format for tximport.
#'
#' @examples
#' object <- S4Vectors::DataFrame(
#'     "txId" = c(
#'         "tx0001",
#'         "tx0002",
#'         "tx0003",
#'         "tx0004"
#'     ),
#'     "geneId" = c(
#'         "gene0001",
#'         "gene0001",
#'         "gene0002",
#'         "gene0002"
#'     )
#' )
#' object <- Tx2Gene(object)
#' con <- file.path(AcidBase::tempdir2(), "tx2gene.csv")
#' export(object = object, con = con)
#' x <- readLines(con, n = 4L)
#' print(x)
#' AcidBase::unlink2(con)
NULL



## Updated 2022-09-13.
`export,Tx2Gene` <- # nolint
    function(object,
             con,
             format, # missing
             ...) {
        if (missing(format)) {
            format <- NULL
        }
        assert(
            validObject(object),
            is.null(format)
        )
        alertInfo(sprintf(
            "Exporting {.cls %s} intentionally without dimnames.",
            "Tx2Gene"
        ))
        export(
            object = as(object, "DFrame"),
            con = con,
            rownames = FALSE,
            colnames = FALSE,
            ...
        )
    }



#' @rdname export
#' @export
setMethod(
    f = "export",
    signature = signature(
        object = "Tx2Gene",
        con = "character",
        format = "missing"
    ),
    definition = `export,Tx2Gene`
)
