#' Import transcript-to-gene annotations
#'
#' Generates a `Tx2Gene` object containing `txId` and `geneId` columns.
#'
#' @note File should not contain column header names.
#' @note Updated 2021-02-01.
#' @export
#'
#' @inheritParams AcidRoxygen::params
#'
#' @param ignoreVersion `logical(2)`.
#'   Ignore transcript ("tx") and/or gene ("gene") versions.
#'
#' @return `Tx2Gene`.
#'
#' @seealso
#' - `stripTranscriptVersions`, `stripGeneVersions`
#'
#' @examples
#' file <- file.path(AcidGenomesTestsURL, "tx2gene.csv")
#' x <- importTx2Gene(
#'     file = file,
#'     organism = "Homo sapiens",
#'     genomeBuild = "GRCh38",
#'     ensemblRelease = 100L
#' )
#' print(x)
importTx2Gene <- function(
    file,
    organism = NULL,
    genomeBuild = NULL,
    release = NULL,
    ignoreVersion = c("tx" = FALSE, "gene" = FALSE)
) {
    assert(
        is.logical(ignoreVersion),
        areSetEqual(
            x = c("tx", "gene"),
            y = names(ignoreVersion)
        )
    )
    data <- import(file = file, rownames = FALSE, colnames = FALSE)
    colnames(data) <- c("txId", "geneId")
    data <- as(data, "DataFrame")
    if (isTRUE(ignoreVersion[["tx"]])) {
        data[["txId"]] <-
            stripTranscriptVersions(data[["txId"]])
    }
    if (isTRUE(ignoreVersion[["gene"]])) {
        data[["geneId"]] <-
            stripGeneVersions(data[["geneId"]])
    }
    metadata(data) <- list(
        "genomeBuild" = genomeBuild,
        "ignoreVersion" = ignoreVersion,
        "organism" = organism,
        "release" = release
    )
    Tx2Gene(data)
}
