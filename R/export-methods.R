## FIXME We need to provide non-breaking support for "file" argument here.
## FIXME pipette colnames = FALSE may not be working the way we want here...



#' @name export
#' @inherit pipette::export description params return title
#' @note Updated 2021-10-21.
#'
#' @details
#' The `Tx2Gene` method automatically disables writing of column names, which
#' is the intended input format for tximport.
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
#' con <- file.path(tempdir(), "tx2gene.csv")
#' export(object = object, con = con)
#' x <- readLines(con, n = 4L)
#' print(x)
#' file.remove(con)
NULL



## Updated 2021-10-21.
`export,Tx2Gene` <-  # nolint
    function(
        object,
        con,
        format,
        ...
    ) {
        stop("FIXME HELLO THERE")
        if (missing(con)) {
            con <- NULL
        }
        if (missing(format)) {
            format <- NULL
        }
        df <- as(object, "DataFrame")
        rownames(df) <- NULL
        export(
            object = df,
            con = con,
            format = format,
            ## FIXME This setting isn't propagating the way we want downstream...
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
