## FIXME Match using HGNC internally for Homo sapiens.
## FIXME Add code coverage against Mus musculus here.



#' Map gene names to Ensembl
#'
#' @export
#' @note Updated 2023-03-03.
#'
#' @details Internally matches using `mapGenesToNCBI` first, so we can support
#' gene synonym matching.
#'
#' @inheritParams AcidRoxygen::params
#'
#' @param genes
#' Gene names (e.g. `"TUT4"`).
#'
#' @examples
#' ## Homo sapiens.
#' x <- mapGeneNamesToEnsembl(
#'     genes = c("TUT4", "ZCCHC11", "TENT3A"),
#'     organism = "Homo sapiens",
#'     genomeBuild = "GRCh38",
#'     release = 109L
#' )
#'
#' ## Mus musculus
#' x <- mapGeneNamesToEnsembl(
#'     genes = c("Nfe2l2", "Nrf2"),
#'     organism = "Mus musculus",
#'     genomeBuild = "GRCm39",
#'     release = 109L
#' )
mapGeneNamesToEnsembl <- function(
        genes,
        organism,
        genomeBuild = NULL,
        release = NULL
    ) {
    assert(
        isCharacter(genes),
        isOrganism(organism),
        isString(genomeBuild, nullOK = TRUE),
        isInt(release, nullOK = TRUE)
    )
    if (is.null(genomeBuild)) {
        genomeBuild <- currentEnsemblGenomeBuild(organism = organism)
    }
    if (is.null(release)) {
        release <- currentEnsemblVersion()
    }
    if (identical(organism, "Homo sapiens")) {
        hgnc <- HGNC()
        map <- as(hgnc, "DataFrame")
        map <- map[, c("hgncId", "ensemblGeneId")]
        ids <- mapGeneNamesToHGNC(genes = genes, hgnc = hgnc)
    } else {
        map <- .importEnsemblNcbiMap(
            organism = organism,
            genomeBuild = genomeBuild,
            release = release
        )
        ids <- mapGeneNamesToNCBI(genes = genes, organism = organism)
    }
    idx <- match(x = ids, table = map[[1L]])
    assert(!anyNA(idx), msg = "Failed to map all genes.")
    out <- map[idx, 2L, drop = TRUE]
    out
}



#' Import Ensembl-to-NCBI gene identifier mappings
#'
#' @note Updated 2023-03-02.
#' @noRd
.importEnsemblNcbiMap <- function(organism, genomeBuild, release) {
    ## Ensure we remove the patch version.
    genomeBuild <- sub(
        pattern = "\\.p[0-9]+$",
        replacement = "",
        x = genomeBuild
    )
    url <- pasteURL(
        "ftp.ensembl.org",
        "pub",
        paste0("release-", release),
        "tsv",
        snakeCase(organism),
        paste0(
            gsub(
                pattern = " ",
                replacement = "_",
                x = organism,
                fixed = TRUE
            ),
            ".",
            genomeBuild,
            ".",
            release,
            ".entrez.tsv.gz"
        ),
        protocol = "https"
    )
    df <- import(con = .cacheIt(url))
    df <- as(df, "DataFrame")
    rownames(df) <- NULL
    colnames(df) <- camelCase(colnames(df))
    keep <- df[["dbName"]] == "EntrezGene"
    df <- df[keep, , drop = FALSE]
    keep <- df[["infoType"]] == "DEPENDENT"
    df <- df[keep, , drop = FALSE]
    df <- df[, c("xref", "geneStableId")]
    df <- df[complete.cases(df), , drop = FALSE]
    df <- unique(df)
    assert(all(grepl(pattern = "^[0-9]+$", x = df[["xref"]])))
    colnames(df) <- c("ncbiGeneId", "ensemblGeneId")
    df[["ncbiGeneId"]] <- as.integer(df[["ncbiGeneId"]])
    idx <- order(df[["ncbiGeneId"]], df[["ensemblGeneId"]])
    df <- df[idx, , drop = FALSE]
    keep <- !duplicated(df[["ncbiGeneId"]])
    df <- df[keep, , drop = FALSE]
    df
}