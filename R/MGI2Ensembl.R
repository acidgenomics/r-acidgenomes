#' @inherit MGI2Ensembl-class title description return
#' @note Updated 2021-02-01.
#' @export
#' @examples
#' object <- MGI2Ensembl()
#' print(object)
MGI2Ensembl <- function() {  # nolint
    assert(hasInternet())
    alert("Importing MGI-to-Ensembl gene ID mappings.")
    url <- pasteURL(
        "www.informatics.jax.org",
        "downloads",
        "reports",
        "MGI_Gene_Model_Coord.rpt",
        protocol = "http"
    )
    file <- .cacheIt(url)
    ## NOTE vroom now warns about expected parsing issues, which cannot be
    ## suppressed.
    ## See related issue:
    ## https://github.com/r-lib/vroom/issues/300
    df <- import(file = file, format = "tsv", colnames = TRUE)
    df <- as(df[, c(1L, 11L)], "DataFrame")
    colnames(df) <- c("mgi", "ensembl")
    df <- df[complete.cases(df), , drop = FALSE]
    df[["mgi"]] <- as.integer(gsub("^MGI\\:", "", df[["mgi"]]))
    assert(hasNoDuplicates(df[["mgi"]]))
    rownames(df) <- df[["mgi"]]
    df <- df[order(df[["mgi"]]), , drop = FALSE]
    new(Class = "MGI2Ensembl", df)
}
