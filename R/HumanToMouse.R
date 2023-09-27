#' Import human-to-mouse gene mappings
#'
#' @export
#' @note Updated 2023-09-26.
#'
#' @param unique `logical(1)`.
#' Only return unique 1:1 mappings between human and mouse.
#' Recommended by default.
#'
#' @return `HumanToMouse`.
#'
#' @seealso
#' - https://www.informatics.jax.org/
#' - https://www.biostars.org/p/9567892/
#'
#' @examples
#' object <- HumanToMouse()
#' print(object)
HumanToMouse <- # nolint
    function(unique = TRUE) {
        assert(
            hasInternet(),
            isFlag(unique)
        )
        url <- pasteUrl(
            "www.informatics.jax.org",
            "downloads",
            "reports",
            "HOM_MouseHumanSequence.rpt",
            protocol = "https"
        )
        file <- .cacheIt(url)
        df <- import(file, format = "tsv")
        df <- as(df, "DFrame")
        colnames(df) <- camelCase(colnames(df))
        cols <- c(
            "dbClassKey", "entrezGeneId", "hgncId", "ncbiTaxonId",
            "mouseMgiId", "omimGeneId", "symbol"
        )
        assert(isSubset(cols, colnames(df)))
        df <- df[, cols, drop = FALSE]
        df[["hgncId"]] <- sub(
            pattern = "^HGNC:",
            replacement = "",
            x = df[["hgncId"]]
        )
        df[["hgncId"]] <- as.integer(df[["hgncId"]])
        df[["mouseMgiId"]] <- sub(
            pattern = "^MGI:",
            replacement = "",
            x = df[["mouseMgiId"]]
        )
        df[["mouseMgiId"]] <- as.integer(df[["mouseMgiId"]])
        df[["omimGeneId"]] <- sub(
            pattern = "^OMIM:",
            replacement = "",
            x = df[["omimGeneId"]]
        )
        df[["omimGeneId"]] <- as.integer(df[["omimGeneId"]])
        spl <- split(df, f = df[["ncbiTaxonId"]])
        assert(is(spl, "SplitDFrameList"))
        hs <- spl[["9606"]]
        assert(
            !anyNA(hs[["humanGeneName"]]),
            !anyNA(hs[["humanHgncId"]]),
            !anyNA(hs[["humanNcbiGeneId"]])
        )
        hs[["ncbiTaxonId"]] <- NULL
        hs[["mouseMgiId"]] <- NULL
        hs <- unique(hs)
        colnames(hs)[colnames(hs) == "entrezGeneId"] <- "humanNcbiGeneId"
        colnames(hs)[colnames(hs) == "hgncId"] <- "humanHgncId"
        colnames(hs)[colnames(hs) == "omimGeneId"] <- "humanOmimGeneId"
        colnames(hs)[colnames(hs) == "symbol"] <- "humanGeneName"
        mm <- spl[["10090"]]
        assert(
            !anyNA(mm[["mouseGeneName"]]),
            !anyNA(mm[["mouseMgiId"]]),
            !anyNA(mm[["mouseNcbiGeneId"]])
        )
        mm[["hgncId"]] <- NULL
        mm[["ncbiTaxonId"]] <- NULL
        mm[["omimGeneId"]] <- NULL
        mm <- unique(mm)
        assert(hasNoDuplicates(mm[["dbClassKey"]]))
        colnames(mm)[colnames(mm) == "entrezGeneId"] <- "mouseNcbiGeneId"
        colnames(mm)[colnames(mm) == "symbol"] <- "mouseGeneName"
        df <- leftJoin(x = hs, y = mm, by = "dbClassKey")
        i <- complete.cases(df[, c("humanGeneName", "mouseGeneName")])
        df <- df[i, , drop = FALSE]
        if (isTRUE(unique)) {
            i <- !isDuplicate(df[["humanGeneName"]]) &
                !isDuplicate(df[["mouseGeneName"]])
            df <- df[i, , drop = FALSE]
            assert(
                hasNoDuplicates(df[["humanGeneName"]]),
                hasNoDuplicates(df[["humanHgncId"]]),
                hasNoDuplicates(df[["humanNcbiGeneId"]]),
                hasNoDuplicates(na.omit(df[["humanOmimGeneId"]])),
                hasNoDuplicates(df[["mouseGeneName"]]),
                hasNoDuplicates(df[["mouseMgiId"]]),
                hasNoDuplicates(df[["mouseNcbiGeneId"]])
            )
        }
        i <- order(df[["humanGeneName"]], df[["mouseGeneName"]])
        j <- sort(colnames(df))
        df <- df[i, j, drop = FALSE]
        metadata(df) <- list(
            "date" = Sys.Date(),
            "humanDupes" = sort(dupes(hs[["humanGeneName"]])),
            "mouseDupes" = sort(dupes(mm[["mouseGeneName"]])),
            "packageVersion" = .pkgVersion,
            "unique" = unique,
            "url" = url
        )
        new(Class = "HumanToMouse", df)
    }
