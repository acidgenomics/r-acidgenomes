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
    ## NOTE vroom now messages a warning about expected parsing issues, which
    ## cannot be suppressed. There's some C++ pointer voodoo going on here.
    ## Falling back to using data.table instead, if possible.
    ## See related issue:
    ## https://github.com/r-lib/vroom/issues/300
    engine <- getOption("acid.import.engine")
    if (isInstalled("data.table")) {
        options("acid.import.engine" = "data.table")
    }
    df <- import(file = file, format = "tsv", colnames = TRUE)
    options("acid.import.engine" = engine)
    df <- as(df[, c(1L, 11L)], "DataFrame")
    colnames(df) <- c("mgiId", "ensemblId")
    df <- df[complete.cases(df), , drop = FALSE]
    df[["mgiId"]] <- as.integer(gsub("^MGI\\:", "", df[["mgiId"]]))
    assert(hasNoDuplicates(df[["mgiId"]]))
    rownames(df) <- df[["mgiId"]]
    df <- df[order(df[["mgiId"]]), , drop = FALSE]
    new(Class = "MGI2Ensembl", df)
}
