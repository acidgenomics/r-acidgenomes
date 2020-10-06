#' Import transcript-to-gene annotations
#'
#' Generates a `Tx2Gene` object containing `transcriptID` and `geneID` columns.
#'
#' @note File should not contain column header names.
#' @note Updated 2020-05-10.
#' @export
#'
#' @inheritParams AcidRoxygen::params
#'
#' @param ignoreTxVersion `logical(1)`.
#'   Don't include the transcript version in the identifier.
#'   **Not recommended** by default when handing off salmon or kallisto output
#'   to tximport-DESeq2 workflow.
#' @param ignoreGeneVersion `logical(1)`.
#'   Don't include the gene version in the identifier.
#'   **Recommended** by default when handing off salmon or kallisto output to
#'   tximport-DESeq2 workflow.
#'
#' @return `Tx2Gene`.
#'
#' @seealso
#' - `stripTranscriptVersions`, `stripGeneVersions`
#' - `download-ensembl-genome` script (which generates `tx2gene.csv` output)
#'   defined in koopa shell bootloader.
#'
#' @examples
#' file <- file.path(AcidGenomesTestsURL, "tx2gene.csv")
#' x <- importTx2Gene(
#'     file = file,
#'     organism = "Homo sapiens",
#'     genomeBuild = "GRCh38",
#'     ensemblRelease = 90L
#' )
#' print(x)
importTx2Gene <- function(
    file,
    organism = NULL,
    genomeBuild = NULL,
    ensemblRelease = NULL,
    ignoreTxVersion = FALSE,
    ignoreGeneVersion = TRUE
) {
    assert(
        isFlag(ignoreTxVersion),
        isFlag(ignoreGeneVersion)
    )
    data <- import(file = file, rownames = FALSE, colnames = FALSE)
    colnames(data) <- c("transcriptID", "geneID")
    data <- as(data, "DataFrame")
    if (isTRUE(ignoreTxVersion)) {
        data[["transcriptID"]] <-
            stripTranscriptVersions(data[["transcriptID"]])
    }
    if (isTRUE(ignoreGeneVersion)) {
        data[["geneID"]] <-
            stripGeneVersions(data[["geneID"]])
    }
    metadata(data) <- list(
        organism = as.character(organism),
        genomeBuild = as.character(genomeBuild),
        ensemblRelease = as.integer(ensemblRelease)
    )
    Tx2Gene(data, metadata = TRUE)
}
