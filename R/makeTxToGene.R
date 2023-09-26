#' Make a TxToGene object
#'
#' @section GFF/GTF file:
#'
#' Remote URLs and compressed files are supported.
#'
#' @name makeTxToGene
#' @note Updated 2023-07-31.
#'
#' @inheritParams AcidRoxygen::params
#' @inheritParams params
#'
#' @return `TxToGene`.
#'
#' @examples
#' ## makeTxToGeneFromEnsembl ====
#' x <- makeTxToGeneFromEnsembl(organism = "Homo sapiens")
#' print(x)
#'
#' ## makeTxToGeneFromEnsDb ====
#' if (goalie::isInstalled("EnsDb.Hsapiens.v75")) {
#'     x <- makeTxToGeneFromEnsDb(object = "EnsDb.Hsapiens.v75")
#'     print(x)
#' }
#'
#' ## makeTxToGeneFromGFF ====
#' file <- AcidBase::pasteURL(
#'     "ftp.ensembl.org",
#'     "pub",
#'     "release-102",
#'     "gtf",
#'     "homo_sapiens",
#'     "Homo_sapiens.GRCh38.102.gtf.gz",
#'     protocol = "https"
#' )
#' x <- makeTxToGeneFromGFF(file = file)
#' print(x)
NULL



#' @rdname makeTxToGene
#' @export
## Updated 2023-07-31.
makeTxToGeneFromEnsembl <-
    function(organism,
             genomeBuild = NULL,
             release = NULL,
             ignoreVersion = FALSE) {
        gr <- makeGRangesFromEnsembl(
            organism = organism,
            level = "transcripts",
            genomeBuild = genomeBuild,
            release = release,
            ignoreVersion = ignoreVersion,
            extraMcols = FALSE
        )
        TxToGene(gr)
    }



#' @rdname makeTxToGene
#' @export
## Updated 2023-07-31.
makeTxToGeneFromEnsDb <-
    function(object,
             ignoreVersion = FALSE) {
        gr <- makeGRangesFromEnsDb(
            object = object,
            level = "transcripts",
            ignoreVersion = ignoreVersion,
            extraMcols = FALSE
        )
        TxToGene(gr)
    }



#' @rdname makeTxToGene
#' @export
## Updated 2023-07-31.
makeTxToGeneFromGFF <-
    function(file,
             ignoreVersion = FALSE) {
        gr <- makeGRangesFromGFF(
            file = file,
            level = "transcripts",
            ignoreVersion = ignoreVersion,
            extraMcols = FALSE
        )
        TxToGene(gr)
    }
