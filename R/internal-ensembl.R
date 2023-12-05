#' Assign extra gene metadata columns (mcols) from Ensembl into GRanges
#'
#' @note Updated 2023-12-05.
#' @noRd
#'
#' @param object `GRanges`.
#'
#' @return `GRanges`.
.addEnsemblFtpMcols <-
    function(object, ignoreVersion) {
        assert(
            is(object, "GRanges"),
            isFlag(ignoreVersion)
        )
        provider <- metadata(object)[["provider"]]
        assert(isSubset(provider, c("Ensembl", "GENCODE")))
        if (
            isFALSE(ignoreVersion) &&
            isSubset("geneIdNoVersion", colnames(mcols(object)))
        ) {
            geneIdCol <- "geneIdNoVersion"
        } else {
            geneIdCol <- "geneId"
        }
        if (isSubset(
            x = c(geneIdCol, "description", "geneSynonyms", "ncbiGeneId"),
            y = colnames(mcols(object))
        )) {
            return(object)
        }
        organism <- metadata(object)[["organism"]]
        genomeBuild <- metadata(object)[["genomeBuild"]]
        if (is.null(genomeBuild)) {
            genomeBuild <- currentEnsemblGenomeBuild(organism)
        }
        genomeBuild2 <- sub(
            pattern = "\\.p[0-9]+$",
            replacement = "",
            x = genomeBuild
        )
        if (isSubset(genomeBuild2, "GRCh37")) {
            return(object)
        }
        release <- metadata(object)[["release"]]
        if (is.null(release)) {
            release <- currentEnsemblVersion()
        }
        if (identical(provider, "GENCODE")) {
            release <- mapGencodeToEnsembl(release)
        }
        alert("Downloading extra gene-level metadata from Ensembl.")
        extraMcols <- .ensemblFtpGeneMetadata(
            organism = organism,
            genomeBuild = genomeBuild,
            release = release
        )
        if (is.null(extraMcols)) {
            return(object)
        }
        if (
            isFALSE(ignoreVersion) &&
            identical(geneIdCol, "geneIdNoVersion")
        ) {
            colnames(extraMcols)[
                colnames(extraMcols) == "geneId"
            ] <- "geneIdNoVersion"
        }
        if (isSubset("description", colnames(mcols(object)))) {
            extraMcols[["description"]] <- NULL
        }
        if (isSubset("geneSynonyms", colnames(mcols(object)))) {
            extraMcols[["geneSynonyms"]] <- NULL
        }
        if (isSubset("ncbiGeneId", colnames(mcols(object)))) {
            extraMcols[["ncbiGeneId"]] <- NULL
        }
        mcols <- leftJoin(x = mcols(object), y = extraMcols, by = geneIdCol)
        mcols(object) <- mcols
        object
    }



#' Get a data frame of extra gene-level metadata from Ensembl FTP server
#'
#' @note Updated 2023-10-12.
#' @noRd
#'
#' @return `DFrame`.
#'
#' @examples
#' df <- .ensemblFtpGeneMetadata(
#'     organism = "Homo sapiens",
#'     genomeBuild = "GRCh38",
#'     release = 109L
#' )
#' df <- .ensemblFtpGeneMetadata(
#'     organism = "Mus musculus",
#'     genomeBuild = "GRCm39",
#'     release = 109L
#' )
.ensemblFtpGeneMetadata <-
    function(organism, genomeBuild, release) {
        assert(
            isOrganism(organism),
            isString(genomeBuild),
            isInt(release)
        )
        genomeBuild <- sub(
            pattern = "\\.p[0-9]+$",
            replacement = "",
            x = genomeBuild
        )
        ## Ensembl is prone to timeouts over FTP, so using HTTPS here instead.
        ftpBaseUrl <- pasteUrl(
            "ftp.ensembl.org",
            "pub",
            paste0("release-", release),
            protocol = "https"
        )
        mysqlSubdir <- getUrlDirList(
            url = pasteUrl(ftpBaseUrl, "mysql"),
            pattern = snakeCase(paste(organism, "core", release))
        )
        assert(isString(mysqlSubdir))
        url <- pasteUrl(ftpBaseUrl, "mysql", mysqlSubdir, "gene.txt.gz")
        gene <- tryCatch(
            expr = {
                import(
                    con = .cacheIt(url),
                    format = "tsv",
                    colnames = FALSE,
                    quote = ""
                )
            },
            error = function(e) {
                alertWarning(sprintf("Failed to import {.url %s}.", url))
                NULL
            }
        )
        if (is.null(gene)) {
            return(NULL)
        }
        url <- pasteUrl(
            ftpBaseUrl, "mysql", mysqlSubdir,
            "external_synonym.txt.gz"
        )
        synonym <- tryCatch(
            expr = {
                import(
                    con = .cacheIt(url),
                    format = "tsv",
                    colnames = FALSE,
                    quote = ""
                )
            },
            error = function(e) {
                alertWarning(sprintf("Failed to import {.url %s}.", url))
                NULL
            }
        )
        if (is.null(synonym)) {
            return(NULL)
        }
        url <- pasteUrl(
            ftpBaseUrl, "tsv", snakeCase(organism),
            paste(
                gsub(
                    pattern = " ",
                    replacement = "_",
                    x = organism
                ),
                genomeBuild,
                release,
                "entrez",
                "tsv.gz",
                sep = "."
            )
        )
        entrez <- tryCatch(
            expr = {
                import(
                    con = .cacheIt(url),
                    format = "tsv",
                    colnames = TRUE,
                    quote = ""
                )
            },
            error = function(e) {
                alertWarning(sprintf("Failed to import {.url %s}.", url))
                NULL
            }
        )
        if (is.null(entrez)) {
            return(NULL)
        }
        colMap <- c(
            "mysqlId" = 8L,
            "geneId" = 13L,
            "description" = 10L
        )
        if (identical(genomeBuild, "GRCh38") && identical(release, 87L)) {
            colMap[["geneId"]] <- 14L
            colMap[["description"]] <- 11L
        }
        df1 <- gene
        df1 <- df1[, colMap]
        colnames(df1) <- names(colMap)
        assert(
            is.integer(df1[["mysqlId"]]),
            is.character(df1[["geneId"]]),
            is.character(df1[["description"]])
        )
        df1 <- df1[complete.cases(df1), , drop = FALSE]
        df1 <- as(df1, "DFrame")
        df2 <- synonym
        df2 <- df2[complete.cases(df2), , drop = FALSE]
        df2 <- unique(df2)
        df2 <- split(x = df2, f = df2[[1L]])
        df2 <- lapply(X = df2, FUN = `[[`, 2L)
        df2 <- as.DataFrame(list(
            "mysqlId" = as.integer(names(df2)),
            "geneSynonyms" = unname(df2)
        ))
        df3 <- entrez
        df3 <- df3[, c("gene_stable_id", "xref"), drop = FALSE]
        df3 <- df3[complete.cases(df3), , drop = FALSE]
        keep <- grepl(pattern = "^[0-9]+$", x = df3[["xref"]])
        df3 <- df3[keep, , drop = FALSE]
        df3[["xref"]] <- as.integer(df3[["xref"]])
        df3 <- unique(df3)
        df3 <- split(x = df3, f = df3[[1L]])
        df3 <- lapply(X = df3, FUN = `[[`, 2L)
        df3 <- as.DataFrame(list(
            "geneId" = names(df3),
            "ncbiGeneId" = unname(df3)
        ))
        out <- leftJoin(x = df1, y = df2, by = "mysqlId")
        out <- leftJoin(x = out, y = df3, by = "geneId")
        out[["mysqlId"]] <- NULL
        out <- out[, sort(colnames(out))]
        out
    }
