#' Show an object
#'
#' @name show
#' @author Michael Steinbaugh
#' @note Updated 2023-04-26.
#'
#' @inheritParams AcidRoxygen::params
#'
#' @return Console output.
NULL



## Updated 2023-04-26.
`show,EnsemblGenes` <- # nolint
    function(object) {
        validObject(object)
        showHeader(object)
    }



#' @rdname show
#' @export
setMethod(
    f = "show",
    signature = signature(object = "EnsemblGenes"),
    definition = `show,EnsemblGenes`
)
