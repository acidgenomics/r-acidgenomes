#' Get extra gene metadata from Ensembl
#'
#' @note Updated 2023-04-12.
#' @noRd
#'
#' @examples
#' df <- .ensemblGeneMetadata(
#'     organism = "Homo sapiens",
#'     genomeBuild = "GRCh38",
#'     release = 109L
#' )
.ensemblGeneMetadata <-
    function(organism, genomeBuild, release) {
        assert(
            isOrganism(organism),
            isString(genomeBuild),
            isInt(release)
        )
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
        colnames(df1) <- c("mysqlId", "geneIdNoVersion", "description")
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
        df3 <- as.DataFrame(list(
            "geneIdNoVersion" = names(df3),
            "ncbiGeneId" = unname(df3)
        ))
        out <- leftJoin(x = df1, y = df2, by = "mysqlId")
        out <- leftJoin(x = out, y = df3, by = "geneIdNoVersion")
        out[["mysqlId"]] <- NULL
        out <- out[, sort(colnames(out))]
        out
    }
