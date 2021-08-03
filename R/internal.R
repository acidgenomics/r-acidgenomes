## FIXME Add file caching option.
## FIXME Call to `.cacheIt` here when `cache = TRUE`.

#' Download multiple genome files in a single call
#'
#' @note Updated 2021-08-03.
#' @noRd
#'
#' @return `character`
#'   Destination files.
.downloadURLs <- function(
    urls,
    outputDir,
    cache = FALSE
) {
    assert(
        allAreURLs(urls),
        isString(outputDir),
        isFlag(cache)
    )
    outputDir <- initDir(outputDir)
    destFiles <- vapply(
        X = urls,
        FUN = function(url) {
            file.path(outputDir, basename(url))
        },
        FUN.VALUE = character(1L),
        USE.NAMES = TRUE
    )
    assert(identical(names(urls), names(destFiles)))
    if (isTRUE(cache)) {
        cacheFiles <- .cacheIt(file = urls)
        file.copy(
            from = cacheFiles,
            to = destFiles,
            overwrite = TRUE,
            recursive = FALSE
        )
    } else {
        mapply(
            url = urls,
            destfile = destFiles,
            FUN = download,
            SIMPLIFY = TRUE,
            USE.NAMES = FALSE
        )
    }
    out <- destFiles
    names(out) <- names(urls)
    assert(allAreFiles(out))
    invisible(out)
}



#' Get an internal function from the package
#'
#' @note Updated 2021-03-03.
#' @noRd
.getFun <- function(x) {
    fun <- get(
        x = x,
        envir = asNamespace(.pkgName),
        inherits = FALSE
    )
    assert(is.function(fun))
    fun
}
