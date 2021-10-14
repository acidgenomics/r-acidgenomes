#' @name export
#' @inherit pipette::export description return title
#' @note Updated 2021-10-14.
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
#' file.remove("tx2gene.csv")
NULL



## Updated 2021-09-24.
`export,Tx2Gene` <-  # nolint
    function(
        object,
        con,
        format,
        ...
    ) {
        if (missing(con)) {
            con <- NULL
        }
        if (missing(format)) {
            format <- NULL
        }
        df <- as(object, "DFrame")
        rownames(df) <- NULL
        export(
            object = df,
            con = con,
            format = format,
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
        con = "ANY",
        format = "ANY"
    ),
    definition = `export,Tx2Gene`
)
