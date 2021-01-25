#' Get an internal function from the package
#'
#' @note Updated 2021-01-25.
#' @noRd
.getFun <- function(x) {
    fun <- get(
        x = x,
        envir = asNamespace(packageName()),
        inherits = FALSE
    )
    assert(is.function(fun))
    fun
}
