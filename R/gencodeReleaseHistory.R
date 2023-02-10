#' GENCODE release history
#'
#' @export
#' @note Updated 2023-02-10.
#'
#' @inheritParams AcidRoxygen::params
#'
#' @return `DataFrame`.
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
        url <- pasteURL(
            "www.gencodegenes.org",
            orgShort,
            "releases.html",
            protocol = "https"
        )
        ## NOTE Don't cache the URL here. This helps avoid the edge case
        ## situation of `currentGencodeVersion` returning a newer version.
        html <- rvest::read_html(url)
        table <- rvest::html_table(html)
        df <- as(table[[1L]], "DataFrame")
        colnames(df) <- camelCase(colnames(df), strict = TRUE)
        df
    }
