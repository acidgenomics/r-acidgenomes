#' @inherit MGI2Ensembl-class title description return
#' @note Updated 2020-10-05.
#' @export
#' @examples
#' object <- MGI2Ensembl()
#' print(object)
MGI2Ensembl <- function() {  # nolint
    assert(hasInternet())
    cli_alert("Importing MGI-to-Ensembl gene ID mappings.")
    url <- pasteURL(
        "www.informatics.jax.org",
        "downloads",
        "reports",
        "MGI_Gene_Model_Coord.rpt",
        protocol = "http"
    )
    file <- cacheURL(url = url, pkg = packageName())
    df <- import(file = file, format = "tsv", colnames = TRUE)
    df <- as(df[, c(1L, 11L)], "DataFrame")
    colnames(df) <- c("mgi", "ensembl")
    df <- df[complete.cases(df), ]
    df[["mgi"]] <- as.integer(gsub("^MGI\\:", "", df[["mgi"]]))
    assert(hasNoDuplicates(df[["mgi"]]))
    rownames(df) <- df[["mgi"]]
    df <- df[order(df[["mgi"]]), , drop = FALSE]
    metadata(df) <- .prototypeMetadata
    new(Class = "MGI2Ensembl", df)
}
