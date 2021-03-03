#' Make a Gene2Symbol object
#'
#' @section GFF/GTF file:
#'
#' Remote URLs and compressed files are supported.
#'
#' @name makeGene2Symbol
#' @note Updated 2021-03-03.
#'
#' @inheritParams Gene2Symbol
#' @inheritParams AcidRoxygen::params
#' @inheritParams params
#'
#' @return `Gene2Symbol`.
#'
#' @examples
#' ## makeGene2SymbolFromEnsembl ====
#' x <- makeGene2SymbolFromEnsembl(organism = "Homo sapiens")
#' print(x)
#'
#' ## makeTx2GeneFromEnsDb ====
#' if ("EnsDb.Hsapiens.v75" %in% rownames(installed.packages())) {
#'     x <- makeGene2SymbolFromEnsDb("EnsDb.Hsapiens.v75")
#'     print(x)
#' }
#'
#' ## makeGene2SymbolFromGFF ====
#' file <- pasteURL(
#'     "ftp.ensembl.org",
#'     "pub",
#'     "release-102",
#'     "gtf",
#'     "homo_sapiens",
#'     "Homo_sapiens.GRCh38.102.gtf.gz",
#'     protocol = "ftp"
#' )
#' x <- makeGene2SymbolFromGFF(file)
#' print(x)
NULL



#' @describeIn makeGene2Symbol Make a `Gene2Symbol` object from Ensembl using
#'   an AnnotationHub lookup.
#' @export
## Updated 2021-02-01.
makeGene2SymbolFromEnsembl <-
    function(
        organism,
        genomeBuild = NULL,
        release = NULL,
        ignoreVersion = FALSE,
        format = c("makeUnique", "unmodified", "1:1")
    ) {
        gr <- makeGRangesFromEnsembl(
            organism = organism,
            genomeBuild = genomeBuild,
            release = release,
            ignoreVersion = ignoreVersion,
            level = "genes"
        )
        Gene2Symbol(object = gr, format = match.arg(format))
    }



#' @describeIn makeGene2Symbol Make a `Gene2Symbol` object from an `EnsDb`
#'   object or annotation package.
#' @export
## Updated 2021-02-01.
makeGene2SymbolFromEnsDb <-
    function(
        object,
        ignoreVersion = FALSE,
        format = c("makeUnique", "unmodified", "1:1")
    ) {
        gr <- makeGRangesFromEnsDb(
            object = object,
            ignoreVersion = ignoreVersion,
            level = "genes"
        )
        Gene2Symbol(object = gr, format = match.arg(format))
    }



#' @describeIn makeGene2Symbol Make a `Gene2Symbol` object from a GFF file.
#' @export
## Updated 2020-02-01.
makeGene2SymbolFromGFF <-
    function(
        file,
        ignoreVersion = FALSE,
        format = c("makeUnique", "unmodified", "1:1")
    ) {
        gr <- makeGRangesFromGFF(
            file = file,
            ignoreVersion = ignoreVersion,
            level = "genes"
        )
        Gene2Symbol(object = gr, format = match.arg(format))
    }
