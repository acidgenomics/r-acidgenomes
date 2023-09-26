#' Import Human Genome Organization (HUGO) Gene Nomenclature Committee (HGNC)
#' metadata
#'
#' @export
#' @note Updated 2023-09-26.
#'
#' @return `HGNC`.
#'
#' @seealso
#' - https://www.genenames.org/
#' - https://www.genenames.org/download/statistics-and-files/
#'
#' @examples
#' object <- HGNC()
#' print(object)
HGNC <- # nolint
    function() {
        alert("Importing HGNC complete set.")
        url <- pasteURL(
            "ftp.ebi.ac.uk",
            "pub",
            "databases",
            "genenames",
            "new",
            "tsv",
            "hgnc_complete_set.txt",
            protocol = "https"
        )
        file <- .cacheIt(url)
        df <- import(con = file, format = "tsv")
        df <- as(df, "DFrame")
        colnames(df) <- camelCase(colnames(df), strict = TRUE)
        assert(
            isSubset(
                x = c(
                    "dateApprovedReserved",
                    "dateModified",
                    "dateNameChanged",
                    "dateSymbolChanged",
                    "entrezId",
                    "hgncId"
                ),
                y = colnames(df)
            ),
            hasNoDuplicates(df[["hgncId"]])
        )
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
        df[["dateApprovedReserved"]] <- as.Date(df[["dateApprovedReserved"]])
        df[["dateModified"]] <- as.Date(df[["dateModified"]])
        df[["dateNameChanged"]] <- as.Date(df[["dateNameChanged"]])
        df[["dateSymbolChanged"]] <- as.Date(df[["dateSymbolChanged"]])
        df[["ncbiGeneId"]] <- as.integer(df[["ncbiGeneId"]])
        metadata(df) <- list(
            "date" = Sys.Date(),
            "packageVersion" = .pkgVersion,
            "url" = url
        )
        new("HGNC", df)
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
