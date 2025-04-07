#' Map gene names to Ensembl
#'
#' @export
#' @note Updated 2025-04-07.
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
#' print(x)
#'
#' ## Mus musculus
#' x <- mapGeneNamesToEnsembl(
#'     genes = c("Nfe2l2", "Nrf2"),
#'     organism = "Mus musculus",
#'     genomeBuild = "GRCm39",
#'     release = 109L
#' )
#' print(x)
mapGeneNamesToEnsembl <-
    function(genes,
             organism,
             genomeBuild = NULL,
             release = NULL) {
        assert(
            isCharacter(genes),
            isOrganism(organism),
            isString(genomeBuild, nullOk = TRUE),
            isInt(release, nullOk = TRUE)
        )
        if (
            identical(organism, "Homo sapiens") &&
            is.null(genomeBuild) && is.null(release)
        ) {
            useHgnc <- TRUE
        } else {
            useHgnc <- FALSE
        }
        if (is.null(genomeBuild)) {
            genomeBuild <- currentEnsemblGenomeBuild(organism = organism)
        }
        if (is.null(release)) {
            release <- currentEnsemblVersion()
        }
        if (isTRUE(useHgnc)) {
            hgnc <- Hgnc()
            map <- as(hgnc, "DFrame")
            map <- map[, c("hgncId", "ensemblGeneId")]
            ids <- mapGeneNamesToHgnc(genes = genes, hgnc = hgnc)
        } else {
            map <- .importEnsemblNcbiMap(
                organism = organism,
                genomeBuild = genomeBuild,
                release = release
            )
            ids <- mapGeneNamesToNcbi(genes = genes, organism = organism)
        }
        idx <- match(x = ids, table = map[[1L]])
        if (anyNA(idx)) {
            abort(sprintf(
                "Mapping failure: %s.",
                toInlineString(genes[is.na(idx)])
            ))
        }
        out <- map[idx, 2L, drop = TRUE]
        out
    }



#' Import Ensembl-to-NCBI gene identifier mappings
#'
#' @note Updated 2023-03-02.
#' @noRd
.importEnsemblNcbiMap <-
    function(organism, genomeBuild, release) {
        ## Ensure we remove the patch version.
        genomeBuild <- sub(
            pattern = "\\.p[0-9]+$",
            replacement = "",
            x = genomeBuild
        )
        url <- pasteUrl(
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
        df <- as(df, "DFrame")
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
