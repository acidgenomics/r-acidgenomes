## FIXME Assert that strand contains "+" and "-" as expected.



#' Import Mouse Genome Informatics (MGI) metadata
#'
#' @export
#' @note Updated 2023-09-19.
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
    lines <- import(con = file, format = "lines")
    cn <- strsplit(lines[[1L]], split = "\t")[[1L]]
    assert(hasLength(cn, n = 15L))
    cn <- sub(pattern = "^[0-9]+\\.\\s", replacement = "", x = cn)
    cn <- camelCase(cn)
    cn[cn == "entrezGeneId"] <- "ncbiGeneId"
    cn <- c(cn, "delete")
    lines <- lines[2L:length(lines)]
    con <- textConnection(lines)
    df <- import(
        con = con,
        format = "tsv",
        colnames = cn,
        naStrings = "NA"
    )
    close(con)
    assert(allAreMatchingFixed(x = df[[1L]], pattern = "MGI:"))
    df <- as(df, "DFrame")
    df[["delete"]] <- NULL
    idCol <- "mgiAccessionId"
    assert(hasNoDuplicates(df[[idCol]]))
    df[[idCol]] <- as.integer(sub(
        pattern = "^MGI\\:",
        replacement = "",
        x = df[[idCol]]
    ))
    rownames(df) <- df[[idCol]]
    df <- df[order(df[[idCol]]), , drop = FALSE]
    new(Class = "MGI", df)
}
