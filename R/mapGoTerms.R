#' Map gene ontology (GO) identifiers to term names
#'
#' @export
#' @note Updated 2023-12-13.
#'
#' @section Alternative approach using GO.db package:
#'
#' This supports lookup of specific keys.
#'
#' \preformatted{
#' library(GO.db)
#' keys <- c("GO:0000001", "GO:0000002")
#' object <- select(
#'     x = GO.db,
#'     keys = keys,
#'     columns = c("GOID", "TERM"),
#'     keytype = "GOID"
#' )
#' }
#'
#' @return `DFrame`.
#' Contains `"id"` and `"name"` columns.
#'
#' @seealso
#' - Bioconductor GO.db package.
#' - https://geneontology.org/docs/download-ontology/
#' - https://www.biostars.org/p/9552810/
#'
#' @examples
#' object <- mapGoTerms()
#' print(object)
mapGoTerms <- function() {
    url <- pasteUrl(
        "purl.obolibrary.org",
        "obo",
        "go",
        "go-basic.obo",
        protocol = "https"
    )
    obo <- import(con = .cacheIt(url), format = "obo")
    df <- as.data.frame(obo)
    assert(isSubset(c("id", "name", "obsolete"), colnames(df)))
    i <- !df[["obsolete"]]
    j <- c("id", "name")
    df <- df[i, j, drop = FALSE]
    assert(hasNoDuplicates(df[["id"]]))
    df <- as(df, "DFrame")
    df <- sort(df)
    rownames(df) <- NULL
    metadata(df) <- list(
        "date" = Sys.Date(),
        "url" = url
    )
    df
}
