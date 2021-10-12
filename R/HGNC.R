#' Import HGNC complete set metadata
#'
#' @export
#' @note Updated 2021-03-19.
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
HGNC <-  # nolint
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
            protocol = "ftp"
        )
        file <- .cacheIt(url)
        suppressWarnings({
            df <- import(file, format = "tsv")
        })
        df <- as(df, "DataFrame")
        colnames(df) <- camelCase(colnames(df), strict = TRUE)
        assert(
            isSubset("hgncId", colnames(df)),
            hasNoDuplicates(df[["hgncId"]])
        )
        df[["hgncId"]] <- as.integer(gsub("^HGNC\\:", "", df[["hgncId"]]))
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
