#' Current genome version
#'
#' Obtain the latest release version from various genome annotation sources.
#'
#' @name currentGenomeVersion
#' @note Updated 2023-09-15.
#'
#' @inheritParams AcidRoxygen::params
#' @param dmel `logical(1)`.
#' Return *Drosophila melanogaster* genome release version.
#'
#' @return `character(1)` or `integer(1)` when possible.
#' Release version.
#'
#' @seealso
#' Refer to the koopa package for Bash shell variants.
#'
#' @examples
#' ## Ensembl.
#' x <- try(
#'     expr = {
#'         currentEnsemblVersion()
#'     },
#'     silent = TRUE
#' )
#' print(x)
#' ## GENCODE.
#' x <- try(
#'     expr = {
#'         currentGencodeVersion(organism = "Homo sapiens")
#'     },
#'     silent = TRUE
#' )
#' print(x)
#' ## RefSeq.
#' x <- try(
#'     expr = {
#'         currentRefseqVersion()
#'     },
#'     silent = TRUE
#' )
#' print(x)
#' ## WormBase.
#' x <- try(
#'     expr = {
#'         currentWormbaseVersion()
#'     },
#'     silent = TRUE
#' )
#' print(x)
NULL



#' @rdname currentGenomeVersion
#' @export
currentEnsemblVersion <- function() {
    x <- import(
        con = pasteUrl(
            "ftp.ensembl.org",
            "pub",
            "current_README",
            protocol = "https"
        ),
        format = "lines",
        quiet = TRUE
    )
    x <- x[[3L]]
    x <- strsplit(x = x, split = " ", fixed = TRUE)[[1L]][[3L]]
    x <- as.integer(x)
    x
}



#' @rdname currentGenomeVersion
#' @export
currentGencodeVersion <-
    function(organism = c("Homo sapiens", "Mus musculus")) {
        assert(requireNamespaces("RCurl"))
        organism <- match.arg(organism)
        url <- "https://www.gencodegenes.org"
        switch(
            EXPR = organism,
            "Homo sapiens" = {
                shortName <- "human"
                pattern <- "Release [[:digit:]]+"
            },
            "Mus musculus" = {
                shortName <- "mouse"
                pattern <- "Release M[[:digit:]]+"
            }
        )
        url <- paste0(url, "/", shortName, "/")
        x <- RCurl::getURL(url)
        x <- strsplit(x = x, split = "\n")[[1L]]
        x <- grep(
            pattern = paste0("^<h1>", pattern),
            x = x,
            value = TRUE
        )
        x <- strsplit(x = x, split = " ", fixed = TRUE)[[1L]][[2L]]
        if (identical(organism, "Homo sapiens")) {
            x <- as.integer(x)
        }
        x
    }



#' @rdname currentGenomeVersion
#' @export
currentRefseqVersion <- function() {
    x <- import(
        con = "https://ftp.ncbi.nlm.nih.gov/refseq/release/RELEASE_NUMBER",
        format = "lines",
        quiet = TRUE
    )
    x <- as.integer(x)
    x
}



#' @rdname currentGenomeVersion
#' @export
currentFlybaseVersion <- function(dmel = FALSE) {
    assert(isFlag(dmel))
    url <- pasteUrl("ftp.flybase.net", "releases", protocol = "ftp")
    if (isTRUE(dmel)) {
        x <- getUrlDirList(pasteUrl(url, "current"))
        x <- grep(pattern = "^dmel_r[.0-9]+$", x = x, value = TRUE)
        x <- strsplit(x = x, split = "_", fixed = TRUE)[[1L]][[2L]]
    } else {
        x <- getUrlDirList(url)
        x <- grep(pattern = "^FB[0-9]{4}_[0-9]{2}$", x = x, value = TRUE)
        x <- tail(sort(x), n = 1L)
    }
    x
}



#' @rdname currentGenomeVersion
#' @export
currentWormbaseVersion <- function() {
    url <- pasteUrl(
        "ftp.wormbase.org", "pub", "wormbase",
        "releases", "current-production-release",
        protocol = "ftp"
    )
    x <- getUrlDirList(paste0(url, "/"))
    x <- grep(pattern = "letter.WS[0-9]+", x = x, value = TRUE)
    x <- strsplit(x = x, split = ".", fixed = TRUE)[[1L]][[2L]]
    x
}
