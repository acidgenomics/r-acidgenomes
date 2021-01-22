#' Get the directives from a GFF file
#'
#' @export
#' @note Updated 2021-01-21.
#'
#' @inheritParams AcidRoxygen::params
#'
#' @details
#' Matches lines beginning with `#!<key> <value>` or `##<key>: <value>`
#'
#' @section GFF3:
#'
#' Lines beginning with '##' are directives (sometimes called pragmas or
#' meta-data) and provide meta-information about the document as a whole. Blank
#' lines should be ignored by parsers and lines beginning with a single '#' are
#' used for human-readable comments and can be ignored by parsers. End-of-line
#' comments (comments preceded by # at the end of and on the same line as a
#' feature or directive line) are not allowed.
#'
#' @seealso
#' - https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
#'
#' @return `DataFrame` or `NULL`.
#'
#' @examples
#' url <- pasteURL(
#'     "ftp.ensembl.org",
#'     "pub",
#'     "release-102",
#'     "gtf",
#'     "homo_sapiens",
#'     "Homo_sapiens.GRCh38.102.gtf.gz",
#'     protocol = "ftp"
#' )
#' df <- getGFFDirectives(url)
#' print(df)
getGFFDirectives <- function(file, nMax = Inf) {
    file <- .cacheIt(file)
    lines <- import(
        file = file,
        format = "lines",
        comment = "",
        nMax = nMax,
        quiet = TRUE
    )
    pattern <- "^(#!|#+)([a-z-]+)(:)?\\s+(.+)$"
    lines <- grep(pattern = pattern, x = lines, value = TRUE)
    if (!hasLength(lines)) return(NULL)
    mat <- str_match(
        string = grep(pattern = pattern, x = lines, value = TRUE),
        pattern = pattern
    )
    assert(is.matrix(mat), hasRows(mat))
    df <- as(mat, "DataFrame")
    df <- df[, c(3L, 5L), drop = FALSE]
    colnames(df) <- c("key", "value")
    df
}
