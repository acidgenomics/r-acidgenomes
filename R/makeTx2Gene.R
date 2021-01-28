#' Make a Tx2Gene object
#'
#' @section GFF/GTF file:
#'
#' Remote URLs and compressed files are supported.
#'
#' @name makeTx2Gene
#' @note Updated 2021-01-28.
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
#' file <- pasteURL(
#'     "ftp.ensembl.org",
#'     "pub",
#'     "release-102",
#'     "gtf",
#'     "homo_sapiens",
#'     "Homo_sapiens.GRCh38.102.gtf.gz",
#'     protocol = "ftp"
#' )
#' x <- makeTx2GeneFromGFF(file)
#' print(x)
NULL



#' @rdname makeTx2Gene
#' @export
## Updated 2021-01-25.
makeTx2GeneFromEnsembl <-
    function() {
        gr <- do.call(
            what = makeGRangesFromEnsembl,
            args = matchArgsToDoCall(
                args = list(
                    "level" = "transcripts",
                    "synonyms" = FALSE
                )
            )
        )
        Tx2Gene(gr)
    }

f <- formals(makeGRangesFromEnsembl)
f <- f[setdiff(names(f), c("level", "synonyms"))]
formals(makeTx2GeneFromEnsembl) <- f



#' @rdname makeTx2Gene
#' @export
## Updated 2021-01-25.
makeTx2GeneFromEnsDb <- function(object) {
    gr <- do.call(
        what = makeGRangesFromEnsDb,
        args = matchArgsToDoCall(
            args = list(
                "level" = "transcripts",
                "synonyms" = FALSE
            )
        )
    )
    Tx2Gene(gr)
}

f <- formals(makeGRangesFromEnsDb)
f <- f[setdiff(names(f), c("level", "synonyms"))]
formals(makeTx2GeneFromEnsDb) <- f



#' @rdname makeTx2Gene
#' @export
## Updated 2021-01-25.
makeTx2GeneFromGFF <- function(file) {
    gr <- do.call(
        what = makeGRangesFromGFF,
        args = matchArgsToDoCall(
            args = list(
                "level" = "transcripts",
                "synonyms" = FALSE
            )
        )
    )
    Tx2Gene(gr)
}

f <- formals(makeGRangesFromGFF)
f <- f[setdiff(names(f), c("level", "synonyms"))]
formals(makeTx2GeneFromGFF) <- f
