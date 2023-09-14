#' NCBI gene history
#'
#' @export
#' @note Updated 2023-09-14.
#'
#' @param `integer(1)`.
#' NCBI taxonomy identifier.
#' For *Homo sapiens*, use `9606L`.
#'
#' @return `DFrame`.
#'
#' @examples
#' ## Homo sapiens.
#' df <- NcbiGeneHistory(taxonomyId = 9606L)
NcbiGeneHistory <- function(taxonomyId) {
    assert(isInt(taxonomyId))
    url <- pasteURL(
        "ftp.ncbi.nih.gov",
        "gene",
        "DATA",
        "gene_history.gz",
        protocol = "ftp"
    )
    df <- import(
        con = .cacheIt(url),
        format = "tsv",
        engine = "readr"
    )
    colnames(df) <- camelCase(colnames(df))
    assert(identical(
        x = colnames(df),
        y = c(
            "xTaxId",
            "geneId",
            "discontinuedGeneId",
            "discontinuedSymbol",
            "discontinueDate"
        )
    ))
    colnames(df)[colnames(df) == "xTaxId"] <- "taxonomyId"
    keep <- df[["taxonomyId"]] == taxonomyId
    df <- df[keep, , drop = FALSE]
    df <- as(df, "DFrame")
    df[["geneId"]] <- as.integer(df[["geneId"]])
    df[["discontinuedGeneId"]] <- as.integer(df[["discontinuedGeneId"]])
    df[["discontinueDate"]] <-
        as.Date(x = as.character(xxx), tryFormats = "%Y%m%d")
    df
}
