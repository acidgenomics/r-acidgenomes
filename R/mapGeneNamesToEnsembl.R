## FIXME Can use dbXrefs metadata in NcbiGeneInfo object for this task.
## FIXME Can we add support for mapping from MGI here (Mus musculus)



#' Map gene names to Ensembl
#'
#' @export
#' @note Updated 2025-04-14.
#'
#' @details Internally matches using `mapGeneNamesToNcbi`, so we can support
#' gene synonym matching.
#'
#' @inheritParams AcidRoxygen::params
#'
#' @param genes
#' Gene names (e.g. `"TUT4"`).
#'
#' @param hgnc `Hgnc` object.
#' Supported for *Homo sapiens* genome only. Snapshot of HGNC annotations.
#' When defined, overrides default behavior of mapping internally with
#' `mapGeneNamesToNcbi` to favor `mapGeneNamesToHgnc` instead.
#'
#' @param ncbi `NcbiGeneInfo` object.
#' Snapshot of NCBI annotations. Passes to `mapGeneNamesToNcbi` internally.
#'
#' @examples
#' ## Homo sapiens.
#' x <- mapGeneNamesToEnsembl(
#'     genes = c("TUT4", "ZCCHC11", "TENT3A"),
#'     organism = "Homo sapiens"
#' )
#' print(x)
#'
#' ## Mus musculus
#' x <- mapGeneNamesToEnsembl(
#'     genes = c("Nfe2l2", "Nrf2"),
#'     organism = "Mus musculus"
#' )
#' print(x)
mapGeneNamesToEnsembl <-
    function(genes,
             organism,
             ignoreCase = FALSE,
             hgnc = NULL,
             ncbi = NULL) {
        assert(
            isCharacter(genes),
            isOrganism(organism),
            isFlag(ignoreCase),
            is(ncbi, "NcbiGeneInfo") || is.null(ncbi),
            is(hgnc, "Hgnc") || is.null(hgnc)
        )
        ## Default to HGNC over NCBI for Homo sapiens.
        if (
            identical(organism, "Homo sapiens") &&
            is.null(hgnc) &&
            is.null(ncbi)
        ) {
            hgnc <- Hgnc()
        }
        if (is(hgnc, "Hgnc")) {
            assert(
                identical(organism, "Homo sapiens"),
                is.null(ncbi)
            )
            map <- as(hgnc, "DFrame")
            map <- map[, c("hgncId", "ensemblGeneId")]
            ids <- mapGeneNamesToHgnc(
                genes = genes,
                ignoreCase = ignoreCase,
                hgnc = hgnc
            )
        } else {
            ## FIXME Rework this to extract from "dbXrefs" column instead.
            stop("FIXME Reworking")
            map <- .importEnsemblNcbiMap(
                organism = organism,
                genomeBuild = genomeBuild,
                release = release
            )
            ids <- mapGeneNamesToNcbi(
                genes = genes,
                organism = organism,
                ignoreCase = ignoreCase,
                ncbi = ncbi
            )
        }
        idx <- match(x = ids, table = map[[1L]])
        if (anyNA(idx)) {
            fail <- genes[is.na(idx)]
            abort(sprintf(
                "%d mapping %s: %s.",
                length(fail),
                ngettext(
                    n = length(fail),
                    msg1 = "failure",
                    msg2 = "failures"
                ),
                toInlineString(x = fail, n = 20L)
            ))
        }
        out <- map[idx, 2L, drop = TRUE]
        out
    }



## FIXME This approach misses some mappings that are defined in NCBI-to-Ensembl
## so need to rethink this approach.

#' Import Ensembl-to-NCBI gene identifier mappings
#'
#' @note Updated 2025-04-15.
#' @noRd
.importEnsemblToNcbiGeneMap <-
    function(
        organism,
        genomeBuild = NULL,
        release = NULL,
        uniqueOnly = TRUE
    ) {
        assert(
            isOrganism(organism),
            isString(genomeBuild, nullOk = TRUE),
            isInt(release, nullOk = TRUE),
            isFlag(uniqueOnly)
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
        df <- df[, c("geneStableId", "xref")]
        df <- df[complete.cases(df), , drop = FALSE]
        df <- unique(df)
        assert(all(grepl(pattern = "^[0-9]+$", x = df[["xref"]])))
        colnames(df)[colnames(df) == "geneStableId"] <- "ensemblGeneId"
        colnames(df)[colnames(df) == "xref"] <- "ncbiGeneId"
        df[["ncbiGeneId"]] <- as.integer(df[["ncbiGeneId"]])
        idx <- order(df)
        df <- df[idx, , drop = FALSE]
        if (isTRUE(uniqueOnly)) {
            for (col in seq_len(ncol(df))) {
                keep <- !duplicated(df[[col]])
                df <- df[keep, , drop = FALSE]
            }
        }
        df
    }
