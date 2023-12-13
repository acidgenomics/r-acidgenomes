#' Import RNAcentral identifier mappings
#'
#' @export
#' @note Updated 2023-12-13.
#'
#' @inheritParams AcidRoxygen::params
#'
#' @param database `character(1)`.
#' Database name. Currently supported: `"Ensembl"`, `"GENCODE"`, `"HGNC"`,
#' `"MGI"`, and `"RefSeq"`. HGNC database is *Homo sapiens* specific, and MGI is
#' *Mus musculus* specific.
#'
#' @return `DFrame`.
#' Contains RNAcentral to specified database RNA transcript identifier mappings.
#' Returns `"rnacentralId"`, `"database"`, `"txId"`, `"taxId"`, `"txBiotype"`,
#' and `"geneId"` columns.
#'
#' @seealso
#' - https://rnacentral.org/
#' - https://ftp.ebi.ac.uk/pub/databases/RNAcentral/
#'
#' @examples
#' ## Get RNAcentral to RefSeq identifier mappings.
#' object <- mapRnacentral(organism = "Homo sapiens", database = "refseq")
#' print(object)
mapRnacentral <- function(organism, database) {
    assert(isOrganism(organism))
    database <- match.arg(
        arg = database,
        choices = c(
            "Ensembl",
            "GENCODE",
            "HGNC",
            "MGI",
            "RefSeq"
        )
    )
    database <- switch(
        EXPR = database,
        "Ensembl" = "ensembl",
        "GENCODE" = "ensembl_gencode",
        "HGNC" = "hgnc",
        "MGI" = "mgi",
        "RefSeq" = "refseq"
    )
    taxId <- .mapOrganismToNcbiTaxId(organism)
    url <- pasteUrl(
        "ftp.ebi.ac.uk",
        "pub",
        "databases",
        "RNAcentral",
        "current_release",
        "id_mapping",
        "database_mappings",
        paste0(database, ".tsv"),
        protocol = "ftp"
    )
    df <- import(
        con = .cacheIt(url),
        colnames = c(
            "rnacentralId",
            "database",
            "txId",
            "taxId",
            "txBiotype",
            "geneId"
        )
    )
    assert(
        isSubset(taxId, unique(df[["taxId"]])),
        msg = sprintf(
            "{.var %s} ({.var %s}) not defined in {.url %s}.",
            organism, taxId, url
        )
    )
    i <- df[["taxId"]] == taxId
    df <- df[i, , drop = FALSE]
    df <- as(df, "DFrame")
    metadata(df) <- list(
        "database" = database,
        "date" = Sys.Date(),
        "organism" = organism,
        "url" = url
    )
    df
}
