#' @inherit MGI2Ensembl-class title description return
#' @note Updated 2022-05-04.
#' @export
#' @examples
#' object <- MGI2Ensembl()
#' print(object)
MGI2Ensembl <- function() { # nolint
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
    ## Base import engine returns NAs introduced by coercion.
    suppressWarnings({
        df <- import(
            con = file,
            format = "tsv",
            colnames = TRUE,
            engine = "readr"
        )
    })
    df <- as(df[, c(1L, 11L)], "DataFrame")
    colnames(df) <- c("mgiId", "ensemblId")
    df <- df[complete.cases(df), , drop = FALSE]
    df[["mgiId"]] <- as.integer(gsub("^MGI\\:", "", df[["mgiId"]]))
    assert(hasNoDuplicates(df[["mgiId"]]))
    rownames(df) <- df[["mgiId"]]
    df <- df[order(df[["mgiId"]]), , drop = FALSE]
    new(Class = "MGI2Ensembl", df)
}
