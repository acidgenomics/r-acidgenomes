#' Map Ensembl release to archive URL.
#'
#' @note Updated 2023-09-27.
#' @export
#'
#' @details
#' Requires the rvest package to be installed.
#'
#' @param release `integer(1)` or `character(1)`.
#' Ensembl release (e.g. `100`).
#'
#' @return `character(1)`.
#' Ensembl release URL.
#'
#' @seealso
#' - `biomaRt::listEnsemblArchives()`.
#'
#' @examples
#' try({
#'     mapEnsemblReleaseToURL(release = 100L)
#' })
mapEnsemblReleaseToURL <- function(release) {
    currentUrl <- pasteURL("useast.ensembl.org", protocol = "https")
    if (is.null(release)) {
        return(currentUrl)
    }
    release <- as.character(release)
    assert(
        requireNamespaces("rvest"),
        isString(release)
    )
    url <- pasteURL(
        "useast.ensembl.org",
        "info",
        "website",
        "archives",
        "index.html",
        protocol = "https"
    )
    assert(isAnExistingURL(url))
    html <- rvest::read_html(url)
    ele <- rvest::html_element(html, css = ".archive-box")
    ele <- rvest::html_element(ele, css = ".spaced")
    txt <- rvest::html_text2(ele)
    spl <- strsplit(txt, split = "\n", fixed = TRUE)[[1L]]
    spl <- strSplit(spl, split = ": ", fixed = TRUE)
    df <- as(spl, "DFrame")
    colnames(df) <- c("name", "date")
    df[["version"]] <- sub(
        pattern = "^Ensembl ",
        replacement = "",
        x = df[["name"]]
    )
    df[["currentRelease"]] <- grepl(pattern = "this site", x = df[["date"]])
    df[["date"]] <- strExtract(df[["date"]], pattern = "[A-Za-z]{3} [0-9]{4}")
    df[["url"]] <- unlist(Map(
        version = df[["version"]],
        date = df[["date"]],
        f = function(version, date) {
            switch(
                EXPR = version,
                "GRCh37" = {
                    pasteURL("grch37.ensembl.org", protocol = "https")
                },
                pasteURL(
                    paste0(
                        sub(
                            pattern = " ",
                            replacement = "",
                            x = tolower(date)
                        ),
                        ".archive.ensembl.org"
                    ),
                    protocol = "https"
                )
            )
        },
        USE.NAMES = FALSE
    ))
    df <- df[, c("name", "date", "url", "version", "currentRelease")]
    assert(isSubset(release, df[["version"]]))
    i <- match(x = release, table = df[["version"]])
    isCurrent <- df[i, "currentRelease"]
    if (isTRUE(isCurrent)) {
        return(currentUrl)
    }
    url <- df[i, "url"]
    assert(isAnExistingURL(url))
    url
}
