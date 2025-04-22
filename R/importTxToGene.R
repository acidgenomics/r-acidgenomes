#' Import transcript-to-gene annotations
#'
#' Generates a `TxToGene` object containing `txId` and `geneId` columns.
#'
#' @note File should not contain column header names.
#' @note Updated 2021-02-01.
#' @export
#'
#' @inheritParams AcidRoxygen::params
#'
#' @param ignoreVersion `logical(2)`.
#' Ignore transcript ("tx") and/or gene ("gene") versions.
#'
#' @return `TxToGene`.
#'
#' @seealso
#' - `stripTranscriptVersions`, `stripGeneVersions`
#'
#' @examples
#' file <- file.path(AcidGenomesTestsUrl, "tx2gene.csv")
#' x <- importTxToGene(
#'     file = file,
#'     organism = "Homo sapiens",
#'     genomeBuild = "GRCh38",
#'     release = 100L
#' )
#' print(x)
importTxToGene <-
    function(
        file,
        organism = NULL,
        genomeBuild = NULL,
        release = NULL,
        ignoreVersion = c(
            "tx" = FALSE,
            "gene" = FALSE
        )
    ) {
        assert(
            is.logical(ignoreVersion),
            areSetEqual(
                x = c("tx", "gene"),
                y = names(ignoreVersion)
            )
        )
        data <- import(con = file, rownames = FALSE, colnames = FALSE)
        colnames(data) <- c("txId", "geneId")
        data <- as(data, "DFrame")
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
        TxToGene(data)
    }
