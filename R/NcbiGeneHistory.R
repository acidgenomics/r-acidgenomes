## FIXME Need to class this.
## FIXME Consider putting discontinuedGeneId first and assigning as rowname.



#' NCBI gene history
#'
#' @export
#' @note Updated 2023-09-14.
#'
#' @inheritParams AcidRoxygen::params
#'
#' @return `DFrame`.
#'
#' @examples
#' ## Homo sapiens.
#' df <- NcbiGeneHistory(organism = "Homo sapiens")
NcbiGeneHistory <- function(organism) {
    assert(isOrganism(organism))
    taxId <- .mapOrganismToNcbiTaxId(organism)
    url <- pasteURL(
        "ftp.ncbi.nih.gov", "gene", "DATA", "gene_history.gz",
        protocol = "ftp"
    )
    df <- import(con = .cacheIt(url), format = "tsv", engine = "readr")
    df <- as(df, "DFrame")
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
    keep <- df[["xTaxId"]] == taxId
    df[["xTaxId"]] <- NULL
    df <- df[keep, , drop = FALSE]
    df[["geneId"]] <- as.integer(df[["geneId"]])
    df[["discontinuedGeneId"]] <- as.integer(df[["discontinuedGeneId"]])
    df[["discontinueDate"]] <- as.Date(
        x = as.character(df[["discontinueDate"]]),
        tryFormats = "%Y%m%d"
    )
    df <- df[order(df), , drop = FALSE]
    df
}
