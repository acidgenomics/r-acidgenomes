#' Import HGNC complete set metadata
#'
#' @export
#' @note Updated 2021-01-06.
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
        file <- cacheURL(url = url, pkg = packageName())
        df <- import(file, format = "tsv")
        df <- as(df, "DataFrame")
        colnames(df) <- camelCase(colnames(df))
        assert(
            isSubset("hgncID", colnames(df)),
            hasNoDuplicates(df[["hgncID"]])
        )
        df[["hgncID"]] <- as.integer(gsub("^HGNC\\:", "", df[["hgncID"]]))
        df <- df[order(df[["hgncID"]]), ]
        rownames(df) <- df[["hgncID"]]
        new("HGNC", df)
    }
