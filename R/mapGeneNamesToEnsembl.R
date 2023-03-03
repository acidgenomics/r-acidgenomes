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
#' x <- mapGeneNamesToEnsembl(
#'     genes = c("TUT4", "ZCCHC11", "TENT3A"),
#'     organism = "Homo sapiens",
#'     genomeBuild = "GRCh38",
#'     release = 108L
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
    map <- import(con = .cacheIt(url))
    colnames(map) <- camelCase(colnames(map))
    keep <- map[["dbName"]] == "EntrezGene"
    map <- map[keep, , drop = FALSE]
    keep <- map[["infoType"]] == "DEPENDENT"
    map <- map[keep, , drop = FALSE]
    map <- map[, c("xref", "geneStableId")]
    map <- map[complete.cases(map), , drop = FALSE]
    map <- unique(map)
    assert(all(grepl(pattern = "^[0-9]+$", x = map[["xref"]])))
    colnames(map) <- c("ncbiGeneId", "ensemblGeneId")
    map[["ncbiGeneId"]] <- as.integer(map[["ncbiGeneId"]])
    idx <- order(map[["ncbiGeneId"]], map[["ensemblGeneId"]])
    map <- map[idx, , drop = FALSE]
    keep <- !duplicated(map[["ncbiGeneId"]])
    map <- map[keep, , drop = FALSE]
    ncbiIds <- mapGeneNamesToNCBI(genes = genes, organism = organism)
    idx <- match(x = ncbiIds, table = map[["ncbiGeneId"]])
    assert(!anyNA(idx), msg = "Failed to map all genes.")
    out <- map[idx, "ensemblGeneId", drop = TRUE]
    out
}
