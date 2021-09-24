## FIXME Need to rework this using BiocIO approach.



#' @name export
#' @inherit pipette::export description return title
#' @note Updated 2021-08-07.
#'
#' @details
#' The `Tx2Gene` method automatically disables writing of column names, which
#' is the intended input format for tximport.
#'
#' @inheritParams AcidRoxygen::params
#'
#' @examples
#' object <- DataFrame(
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
#' export(object, file = "tx2gene.csv")
#' unlink("tx2gene.csv")
NULL



## Updated 2021-08-07.
`export,Tx2Gene` <-  # nolint
    function(object, ...) {
        df <- as(object, "DataFrame")
        rownames(df) <- NULL
        export(object = df, colnames = FALSE, ...)
    }



#' @rdname export
#' @export
setMethod(
    f = "export",
    signature = signature("Tx2Gene"),
    definition = `export,Tx2Gene`
)
