#' Import Mouse Genome Informatics (MGI) metadata from the Jackson Laboratory
#'
#' @export
#' @note Updated 2023-03-01.
#'
#' @return `MGI`.
#'
#' @seealso
#' - https://www.informatics.jax.org/
#'
#' @examples
#' object <- MGI()
#' print(object)
MGI <- function() { # nolint
    alert("Importing MGI metadata.")
    url <- pasteURL(
        "www.informatics.jax.org",
        "downloads",
        "reports",
        "MGI_Gene_Model_Coord.rpt",
        protocol = "https"
    )
    file <- .cacheIt(url)
    ## Base import engine returns NAs introduced by coercion.
    suppressWarnings({
        df <- import(
            con = file,
            format = "tsv",
            colnames = TRUE,
            engine = "readr"
        )
    })
    df <- as(df, "DataFrame")
    cn <- colnames(df)
    cn <- sub(pattern = "^X[0-9]+_", replacement = "", x = cn)
    cn <- camelCase(cn)
    cn[cn == "entrezGeneId"] <- "ncbiGeneId"
    colnames(df) <- cn
    idCol <- "mgiAccessionId"
    df[[idCol]] <- as.integer(sub(
        pattern = "^MGI\\:",
        replacement = "",
        x = df[[idCol]]
    ))
    assert(hasNoDuplicates(df[[idCol]]))
    rownames(df) <- df[[idCol]]
    df <- df[order(df[[idCol]]), , drop = FALSE]
    new(Class = "MGI", df)
}
