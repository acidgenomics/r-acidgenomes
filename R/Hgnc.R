#' Import Human Genome Organization (HUGO) Gene Nomenclature Committee (HGNC)
#' metadata
#'
#' @export
#' @note Updated 2025-03-24.
#'
#' @return `Hgnc`.
#'
#' @seealso
#' - https://www.genenames.org/
#' - https://www.genenames.org/download/
#' - https://www.genenames.org/download/archive/
#' - https://www.genenames.org/download/statistics-and-files/
#'
#' @examples
#' object <- Hgnc()
#' print(object)
Hgnc <- # nolint
    function() {
        alert("Importing HGNC complete set.")
        url <- pasteUrl(
            "storage.googleapis.com",
            "public-download-files",
            "hgnc",
            "tsv",
            "tsv",
            "hgnc_complete_set.txt",
            protocol = "https"
        )
        con <- .cacheIt(url)
        df <- import(con = con, format = "tsv")
        df <- as(df, "DFrame")
        colnames(df) <- camelCase(colnames(df), strict = TRUE)
        assert(
            isSubset(
                x = c(
                    "agr",
                    "dateApprovedReserved",
                    "dateModified",
                    "dateNameChanged",
                    "dateSymbolChanged",
                    "entrezId",
                    "gencc",
                    "hgncId",
                    "name",
                    "symbol"
                ),
                y = colnames(df)
            ),
            hasNoDuplicates(df[["hgncId"]])
        )
        df[["agr"]] <- NULL
        df[["gencc"]] <- NULL
        df[["hgncId"]] <- as.integer(sub(
            pattern = "^HGNC\\:",
            replacement = "",
            x = df[["hgncId"]]
        ))
        df <- df[order(df[["hgncId"]]), , drop = FALSE]
        rownames(df) <- df[["hgncId"]]
        isNested <- bapply(
            X = df,
            FUN = function(x) {
                any(grepl(pattern = "|", x = x, fixed = TRUE))
            }
        )
        if (any(isNested)) {
            vars <- names(isNested)[isNested]
            df <- mutateAt(
                object = df,
                vars = vars,
                fun = .splitToCharacterList,
                split = "|"
            )
        }
        colnames(df)[colnames(df) == "entrezId"] <- "ncbiGeneId"
        colnames(df)[colnames(df) == "name"] <- "description"
        colnames(df)[colnames(df) == "symbol"] <- "geneName"
        df[["dateApprovedReserved"]] <- as.Date(df[["dateApprovedReserved"]])
        df[["dateModified"]] <- as.Date(df[["dateModified"]])
        df[["dateNameChanged"]] <- as.Date(df[["dateNameChanged"]])
        df[["dateSymbolChanged"]] <- as.Date(df[["dateSymbolChanged"]])
        df[["ncbiGeneId"]] <- as.integer(df[["ncbiGeneId"]])
        metadata(df) <- list(
            "date" = Sys.Date(),
            "organism" = "Homo sapiens",
            "packageVersion" = .pkgVersion,
            "url" = url
        )
        new("Hgnc", df)
    }



## Updated 2021-03-19.
.splitToCharacterList <- function(x, split = "|") {
    if (all(is.na(x))) {
        return(x)
    }
    x <- strsplit(x, split = split, fixed = TRUE)
    x <- CharacterList(x)
    ## > x <- sort(unique(x))
    x
}
