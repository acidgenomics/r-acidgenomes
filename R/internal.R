#' Download multiple genome files in a single call
#'
#' @note Updated 2021-01-07.
#' @noRd
#'
#' @return `character`
#'   Destination files.
.downloadURLs <- function(urls, outputDir) {
    assert(
        allAreURLs(urls),
        isString(outputDir)
    )
    outputDir <- initDir(outputDir)
    destfiles <- vapply(
        X = urls,
        FUN = function(url) {
            file.path(outputDir, basename(url))
        },
        FUN.VALUE = character(1L)
    )
    assert(identical(names(urls), names(destfiles)))
    out <- mapply(
        url = urls,
        destfile = destfiles,
        FUN = download,
        SIMPLIFY = TRUE,
        USE.NAMES = FALSE
    )
    names(out) <- names(urls)
    assert(allAreFiles(out))
    invisible(out)
}



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
