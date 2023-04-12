## FIXME This is problematic for Mus musculus GRCm39 109...debug.

#' Get extra gene metadata columns (mcols) from Ensembl
#'
#' @note Updated 2023-04-12.
#' @noRd
#'
#' @param object `GRanges`.
#'
#' @return `GRanges`.
#'
#' @examples
#' df <- .ensemblFtpGeneMcols(
#'     organism = "Homo sapiens",
#'     genomeBuild = "GRCh38",
#'     release = 109L
#' )
.addEnsemblFtpMcols <-
    function(object, ignoreVersion) {
        assert(
            is(object, "GRanges"),
            isFlag(ignoreVersion)
        )
        geneIdCol <- ifelse(
            test = ignoreVersion,
            yes = "geneId",
            no = "geneIdNoVersion"
        )
        if (isSubset(
            x = c(
                geneIdCol,
                "description",
                "geneSynonyms",
                "ncbiGeneId"
            ),
            y = colnames(mcols(object))
        )) {
            return(object)
        }
        genomeBuild <- metadata(object)[["genomeBuild"]]
        organism <- metadata(object)[["organism"]]
        release <- metadata(object)[["release"]]
        ftpBaseUrl <- pasteURL(
            "ftp.ensembl.org",
            "pub",
            paste0("release-", release),
            protocol = "ftp"
        )
        mysqlSubdir <- getURLDirList(
            url = pasteURL(ftpBaseUrl, "mysql"),
            pattern = snakeCase(paste(organism, "core", release))
        )
        assert(isString(mysqlSubdir))
        ## FIXME This step is problematic for Mus musculus.
        ## "ftp://ftp.ensembl.org/pub/release-109/mysql/mus_musculus_core_109_39/gene.txt.gz"
        gene <- import(
            con = .cacheIt(pasteURL(
                ftpBaseUrl, "mysql", mysqlSubdir,
                "gene.txt.gz"
            )),
            format = "tsv",
            colnames = FALSE
        )
        synonym <- import(
            con = .cacheIt(pasteURL(
                ftpBaseUrl, "mysql", mysqlSubdir,
                "external_synonym.txt.gz"
            )),
            format = "tsv",
            colnames = FALSE
        )
        entrez <- import(
            con = .cacheIt(pasteURL(
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
            )),
            format = "tsv",
            colnames = TRUE
        )
        df1 <- gene
        df1 <- df1[, c(8, 13, 10)]
        colnames(df1) <- c("mysqlId", geneIdCol, "description")
        df1 <- as(df1, "DataFrame")
        df2 <- synonym
        df2 <- unique(df2)
        df2 <- split(x = df2, f = df2[[1L]])
        df2 <- lapply(X = df2, FUN = `[[`, 2L)
        df2 <- as.DataFrame(list(
            "mysqlId" = names(df2),
            "geneSynonyms" = unname(df2)
        ))
        df3 <- entrez
        df3 <- df3[, c("gene_stable_id", "xref")]
        keep <- grepl(pattern = "^[0-9]+$", x = df3[["xref"]])
        df3 <- df3[keep, ]
        df3[["xref"]] <- as.integer(df3[["xref"]])
        df3 <- unique(df3)
        df3 <- split(x = df3, f = df3[[1L]])
        df3 <- lapply(X = df3, FUN = `[[`, 2L)
        df3lst <- list()
        df3lst[[geneIdCol]] <- names(df3)
        df3lst[["ncbiGeneId"]] <- unname(df3)
        df3 <- as.DataFrame(df3lst)
        mcols <- leftJoin(x = df1, y = df2, by = "mysqlId")
        mcols <- leftJoin(x = mcols, y = df3, by = geneIdCol)
        mcols[["mysqlId"]] <- NULL
        mcols <- mcols[, sort(colnames(mcols))]
        if (isSubset("description", colnames(mcols(object)))) {
            mcols[["description"]] <- NULL
        }
        if (isSubset("geneSynonyms", colnames(mcols(object)))) {
            mcols[["geneSynonyms"]] <- NULL
        }
        if (isSubset("ncbiGeneId", colnames(mcols(object)))) {
            mcols[["ncbiGeneId"]] <- NULL
        }
        mcols <- leftJoin(x = mcols(object), y = df, by = geneIdCol)
        mcols(object) <- mcols
        object
    }
