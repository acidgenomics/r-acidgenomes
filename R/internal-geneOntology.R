#' Map gene ontology (GO) identifiers to term names
#'
#' @note Updated 2023-12-13.
#' @noRd
#'
#' @seealso
#' - https://geneontology.org/docs/download-ontology/
#' - Bioconductor GO.db package
#' - https://www.biostars.org/p/9552810/
.mapGoTerms <- function() {
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
    df <- as(df, "DFrame")
    assert(hasNoDuplicates(df[["id"]]))
    df <- sort(df)
    rownames(df) <- NULL
    metadata(df) <- list(
        "date" = Sys.Date(),
        "url" = url
    )
    df
}
