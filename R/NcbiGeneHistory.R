#' NCBI gene history
#'
#' @export
#' @note Updated 2023-09-15.
#'
#' @inheritParams AcidRoxygen::params
#'
#' @return `NcbiGeneHistory`.
#'
#' @examples
#' ## Homo sapiens.
#' object <- NcbiGeneHistory(organism = "Homo sapiens")
#' print(object)
NcbiGeneHistory <- function(organism) {
    assert(isOrganism(organism))
    taxId <- .mapOrganismToNcbiTaxId(organism)
    assert(isInt(taxId))
    url <- pasteURL(
        "ftp.ncbi.nih.gov", "gene", "DATA", "gene_history.gz",
        protocol = "https"
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
        ),
        hasNoDuplicates(df[["discontinuedGeneId"]])
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
    j <- c(
        "discontinuedGeneId",
        "discontinuedSymbol",
        "discontinueDate",
        "geneId"
    )
    df <- df[, j, drop = FALSE]
    i <- order(df)
    df <- df[i, , drop = FALSE]
    rownames(df) <- df[["discontinuedGeneId"]]
    ## FIXME Include package name and package version.
    metadata(df) <- list(
        "date" = Sys.Date(),
        "organism" = organism,
        "taxonomyId" = taxId,
        "url" = url
    )
    new(df, Class = "NcbiGeneHistory")
}
