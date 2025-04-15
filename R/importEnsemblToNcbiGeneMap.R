#' Import Ensembl-to-NCBI gene identifier mappings
#'
#' @export
#' @note Updated 2025-04-15.
#'
#' @inheritParams AcidRoxygen::params
#'
#' @param uniqueOnly `logical(1)`.
#' Resolve duplicates and only return unique 1:1 mappings.
#'
#' @return `DFrame`.
#'
#' @examples
#' df <- importEnsemblToNcbiGeneMap(organism = "Homo sapiens")
#' print(df)
importEnsemblToNcbiGeneMap <-
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
                keep <- !isDuplicate(df[[col]])
                df <- df[keep, , drop = FALSE]
            }
        }
        df
    }
