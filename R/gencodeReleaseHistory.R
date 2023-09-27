#' GENCODE release history
#'
#' @export
#' @note Updated 2023-09-26.
#'
#' @details
#' Requires the rvest package to be installed.
#'
#' @inheritParams AcidRoxygen::params
#'
#' @return `DFrame`.
#'
#' @examples
#' object <- gencodeReleaseHistory(organism = "Homo sapiens")
gencodeReleaseHistory <-
    function(organism) {
        assert(requireNamespaces("rvest"))
        organism <- match.arg(
            arg = organism,
            choices = c("Homo sapiens", "Mus musculus")
        )
        orgShort <- switch(
            EXPR = organism,
            "Homo sapiens" = "human",
            "Mus musculus" = "mouse"
        )
        url <- pasteUrl(
            "www.gencodegenes.org",
            orgShort,
            "releases.html",
            protocol = "https"
        )
        ## Don't cache the URL here. This helps avoid the edge case situation of
        ## `currentGencodeVersion` returning a newer version.
        html <- rvest::read_html(url)
        table <- rvest::html_table(html)
        df <- as(table[[1L]], "DFrame")
        colnames(df) <- camelCase(colnames(df), strict = TRUE)
        df
    }
