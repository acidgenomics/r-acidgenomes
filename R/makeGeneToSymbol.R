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
#' ## makeGeneToSymbolFromEnsDb ====
#' ## > if (goalie::isInstalled("EnsDb.Hsapiens.v75")) {
#' ## >     x <- makeGeneToSymbolFromEnsDb("EnsDb.Hsapiens.v75")
#' ## >     print(x)
#' ## > }
#'
#' ## makeGeneToSymbolFromGff ====
#' ## > file <- AcidBase::pasteUrl(
#' ## >     "ftp.ensembl.org",
#' ## >     "pub",
#' ## >     "release-102",
#' ## >     "gtf",
#' ## >     "homo_sapiens",
#' ## >     "Homo_sapiens.GRCh38.102.gtf.gz",
#' ## >     protocol = "ftp"
#' ## > )
#' ## > x <- makeGeneToSymbolFromGff(
#' ## >     file = file,
#' ## >     ignoreVersion = FALSE
#' ## > )
#' ## > print(x)
NULL



#' @describeIn makeGeneToSymbol Make a `GeneToSymbol` object from Ensembl using
#' an AnnotationHub lookup.
#' @export
## Updated 2023-12-04.
makeGeneToSymbolFromEnsembl <-
    function(organism,
             genomeBuild = NULL,
             release = NULL,
             ignoreVersion = FALSE,
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
## Updated 2023-12-04.
makeGeneToSymbolFromEnsDb <-
    function(object,
             ignoreVersion = FALSE,
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
## Updated 2023-12-04.
makeGeneToSymbolFromGff <-
    function(file,
             ignoreVersion = FALSE,
             format = c("makeUnique", "1:1", "unmodified")) {
        gr <- makeGRangesFromGff(
            file = file,
            ignoreVersion = ignoreVersion,
            level = "genes"
        )
        GeneToSymbol(object = gr, format = match.arg(format))
    }
