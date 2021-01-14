#' Dynamically cache a URL into package or return local file
#'
#' @note Updated 2021-01-14.
#' @noRd
.cache <- function(file) {
    assert(isString(file))
    if (isAURL(file)) {
        x <- cacheURL(url = file, pkg = packageName())
    } else {
        x <- file
    }
    assert(isAFile(tmpfile))
    x
}
