#' Make a Tx2Gene object
#'
#' @section GFF/GTF file:
#'
#' Remote URLs and compressed files are supported.
#'
#' @name makeTx2Gene
#' @note Updated 2021-03-10.
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
#' if ("EnsDb.Hsapiens.v75" %in% rownames(installed.packages())) {
#'     x <- makeTx2GeneFromEnsDb("EnsDb.Hsapiens.v75")
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
#' x <- makeTx2GeneFromGFF(file = file, ignoreVersion = FALSE)
#' print(x)
NULL



#' @rdname makeTx2Gene
#' @export
## Updated 2021-03-10.
makeTx2GeneFromEnsembl <-
    function(
        organism,
        genomeBuild = NULL,
        release = NULL,
        ignoreVersion = TRUE
    ) {
        gr <- makeGRangesFromEnsembl(
            organism = organism,
            genomeBuild = genomeBuild,
            release = release,
            ignoreVersion = ignoreVersion,
            level = "transcripts"
        )
        Tx2Gene(gr)
    }



#' @rdname makeTx2Gene
#' @export
## Updated 2021-03-10.
makeTx2GeneFromEnsDb <-
    function(
        object,
        ignoreVersion = TRUE
    ) {
        gr <- makeGRangesFromEnsDb(
            object = object,
            ignoreVersion = ignoreVersion,
            level = "transcripts"
        )
        Tx2Gene(gr)
    }



#' @rdname makeTx2Gene
#' @export
## Updated 2021-03-10.
makeTx2GeneFromGFF <-
    function(
        file,
        ignoreVersion = TRUE
    ) {
        gr <- makeGRangesFromGFF(
            file = file,
            ignoreVersion = ignoreVersion,
            level = "transcripts"
        )
        Tx2Gene(gr)
    }
