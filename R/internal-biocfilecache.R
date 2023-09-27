#' Dynamically cache a URL into package or return local file
#'
#' @note Updated 2021-08-03.
#' @note R won't let you name the function `.cache`.
#' @noRd
#'
#' @param file `character`.
#' URLs and/or local files.
#'
#' @return `character`.
#' Local file paths.
#'
#' @examples
#' file <- pasteUrl(AcidGenomesTestsUrl, "ensembl.gtf")
#' tmpfile <- .cacheIt(file)
#' print(tmpfile)
.cacheIt <- function(file) {
    vapply(
        X = file,
        FUN = function(file) {
            if (isAUrl(file)) {
                x <- cacheUrl(url = file, pkg = .pkgName)
            } else {
                x <- file
            }
            assert(isAFile(x))
            x
        },
        FUN.VALUE = character(1L),
        USE.NAMES = FALSE
    )
}
