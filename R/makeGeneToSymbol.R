#' Make a GeneToSymbol object
#'
#' @section GFF/GTF file:
#'
#' Remote URLs and compressed files are supported.
#'
#' @name makeGeneToSymbol
#' @note Updated 2021-08-03.
#'
#' @inheritParams GeneToSymbol
#' @inheritParams AcidRoxygen::params
#' @inheritParams params
#'
#' @return `GeneToSymbol`.
#'
#' @examples
#' ## makeGeneToSymbolFromEnsembl ====
#' x <- makeGeneToSymbolFromEnsembl(
#'     organism = "Homo sapiens",
#'     ignoreVersion = FALSE
#' )
#' print(x)
#'
#' ## makeTxToGeneFromEnsDb ====
#' if (goalie::isInstalled("EnsDb.Hsapiens.v75")) {
#'     x <- makeGeneToSymbolFromEnsDb("EnsDb.Hsapiens.v75")
#'     print(x)
#' }
#'
#' ## makeGeneToSymbolFromGFF ====
#' file <- AcidBase::pasteURL(
#'     "ftp.ensembl.org",
#'     "pub",
#'     "release-102",
#'     "gtf",
#'     "homo_sapiens",
#'     "Homo_sapiens.GRCh38.102.gtf.gz",
#'     protocol = "ftp"
#' )
#' x <- makeGeneToSymbolFromGFF(
#'     file = file,
#'     ignoreVersion = FALSE
#' )
#' print(x)
NULL



#' @describeIn makeGeneToSymbol Make a `GeneToSymbol` object from Ensembl using
#' an AnnotationHub lookup.
#' @export
## Updated 2021-08-03.
makeGeneToSymbolFromEnsembl <-
    function(organism,
             genomeBuild = NULL,
             release = NULL,
             ignoreVersion = TRUE,
             format = c("makeUnique", "1:1", "unmodified")) {
        gr <- makeGRangesFromEnsembl(
            organism = organism,
            genomeBuild = genomeBuild,
            release = release,
            ignoreVersion = ignoreVersion,
            level = "genes"
        )
        GeneToSymbol(object = gr, format = match.arg(format))
    }



#' @describeIn makeGeneToSymbol Make a `GeneToSymbol` object from an `EnsDb`
#' object or annotation package.
#' @export
## Updated 2021-08-03.
makeGeneToSymbolFromEnsDb <-
    function(object,
             ignoreVersion = TRUE,
             format = c("makeUnique", "1:1", "unmodified")) {
        gr <- makeGRangesFromEnsDb(
            object = object,
            ignoreVersion = ignoreVersion,
            level = "genes"
        )
        GeneToSymbol(object = gr, format = match.arg(format))
    }



#' @describeIn makeGeneToSymbol Make a `GeneToSymbol` object from a GFF file.
#' @export
## Updated 2020-08-03.
makeGeneToSymbolFromGFF <-
    function(file,
             ignoreVersion = TRUE,
             format = c("makeUnique", "1:1", "unmodified")) {
        gr <- makeGRangesFromGFF(
            file = file,
            ignoreVersion = ignoreVersion,
            level = "genes"
        )
        GeneToSymbol(object = gr, format = match.arg(format))
    }
