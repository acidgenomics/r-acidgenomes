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
    assert(identical(
        x = colnames(df),
        y = c(
            "X_tax_id",
            "GeneID",
            "Discontinued_GeneID",
            "Discontinued_Symbol",
            "Discontinue_Date"
        )
    ))
    colnames(df) <- camelCase(colnames(df))
    colnames(df)[colnames(df) == "xTaxId"] <- "taxonomyId"
    keep <- df[["taxonomyId"]] == taxonomyId
    df <- df[keep, , drop = FALSE]
    df <- as(df, "DFrame")
    df
}
