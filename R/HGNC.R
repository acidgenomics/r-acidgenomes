#' Import HGNC complete set metadata
#'
#' @export
#' @note Updated 2021-01-14.
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
        df <- import(file, format = "tsv")
        df <- as(df, "DataFrame")
        colnames(df) <- camelCase(colnames(df), strict = TRUE)
        assert(
            isSubset("hgncId", colnames(df)),
            hasNoDuplicates(df[["hgncId"]])
        )
        df[["hgncId"]] <- as.integer(gsub("^HGNC\\:", "", df[["hgncId"]]))
        df <- df[order(df[["hgncId"]]), ]
        rownames(df) <- df[["hgncId"]]
        new("HGNC", df)
    }
