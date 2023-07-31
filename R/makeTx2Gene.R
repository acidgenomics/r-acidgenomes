#' Make a Tx2Gene object
#'
#' @section GFF/GTF file:
#'
#' Remote URLs and compressed files are supported.
#'
#' @name makeTx2Gene
#' @note Updated 2023-07-31.
#'
#' @inheritParams AcidRoxygen::params
#' @inheritParams params
#'
#' @return `Tx2Gene`.
#'
#' @examples
#' ## makeTx2GeneFromEnsembl ====
#' x <- makeTx2GeneFromEnsembl(organism = "Homo sapiens")
#' print(x)
#'
#' ## makeTx2GeneFromEnsDb ====
#' if (goalie::isInstalled("EnsDb.Hsapiens.v75")) {
#'     x <- makeTx2GeneFromEnsDb(object = "EnsDb.Hsapiens.v75")
#'     print(x)
#' }
#'
#' ## makeTx2GeneFromGFF ====
#' file <- AcidBase::pasteURL(
#'     "ftp.ensembl.org",
#'     "pub",
#'     "release-102",
#'     "gtf",
#'     "homo_sapiens",
#'     "Homo_sapiens.GRCh38.102.gtf.gz",
#'     protocol = "ftp"
#' )
#' x <- makeTx2GeneFromGFF(file = file)
#' print(x)
NULL



#' @rdname makeTx2Gene
#' @export
## Updated 2023-07-31.
makeTx2GeneFromEnsembl <-
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
        Tx2Gene(gr)
    }



#' @rdname makeTx2Gene
#' @export
## Updated 2023-07-31.
makeTx2GeneFromEnsDb <-
    function(object,
             ignoreVersion = FALSE) {
        gr <- makeGRangesFromEnsDb(
            object = object,
            level = "transcripts",
            ignoreVersion = ignoreVersion,
            extraMcols = FALSE
        )
        Tx2Gene(gr)
    }



#' @rdname makeTx2Gene
#' @export
## Updated 2023-07-31.
makeTx2GeneFromGFF <-
    function(file,
             ignoreVersion = FALSE) {
        gr <- makeGRangesFromGFF(
            file = file,
            level = "transcripts",
            ignoreVersion = ignoreVersion,
            extraMcols = FALSE
        )
        Tx2Gene(gr)
    }
