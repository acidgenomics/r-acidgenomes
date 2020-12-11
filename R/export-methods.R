#' @name export
#' @inherit pipette::export description return title
#' @note Updated 2020-12-10.
#'
#' @details
#' The `Tx2Gene` method automatically disables writing of column names, which
#' is the intended input format for tximport.
#'
#' @inheritParams AcidRoxygen::params
#'
#' @examples
#' object <- DataFrame(
#'     transcriptID = c(
#'         "tx0001",
#'         "tx0002",
#'         "tx0003",
#'         "tx0004"
#'     ),
#'     geneID = c(
#'         "gene0001",
#'         "gene0001",
#'         "gene0002",
#'         "gene0002"
#'     )
#' )
#' export(object, file = "tx2gene.csv")
#' unlink("tx2gene.csv")
NULL



## Updated 2020-12-10.
`export,Tx2Gene` <-  # nolint
    function(object, ...) {
        export(
            object = as(object, "DataFrame"),
            colnames = FALSE,
            ...
        )
    }



#' @rdname export
#' @export
setMethod(
    f = "export",
    signature = signature("Tx2Gene"),
    definition = `export,Tx2Gene`
)
