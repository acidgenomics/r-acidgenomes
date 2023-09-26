## FIXME Need to improve formatting of dates.



#' Import Human Genome Organization (HUGO) Gene Nomenclature Committee (HGNC)
#' metadata
#'
#' @export
#' @note Updated 2023-09-15.
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
        ## FIXME Switch to base engine here.
        df <- import(con = file, format = "tsv")
        df <- as(df, "DFrame")
        colnames(df) <- camelCase(colnames(df), strict = TRUE)
        idCol <- "hgncId"
        assert(
            isSubset(idCol, colnames(df)),
            hasNoDuplicates(df[[idCol]])
        )
        df[[idCol]] <- as.integer(sub(
            pattern = "^HGNC\\:",
            replacement = "",
            x = df[[idCol]]
        ))
        df <- df[order(df[[idCol]]), , drop = FALSE]
        rownames(df) <- df[[idCol]]
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
