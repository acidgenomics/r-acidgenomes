#' Make a Tx2Gene object
#'
#' @section GFF/GTF file:
#'
#' Remote URLs and compressed files are supported.
#'
#' @name makeTx2Gene
#' @note Updated 2021-01-08.
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
#' ## GTF
#' file <- file.path(AcidGenomesTestsURL, "example.gtf")
#' x <- makeTx2GeneFromGFF(file)
#' print(x)
#'
#' ## GFF3
#' file <- file.path(AcidGenomesTestsURL, "example.gff3")
#' x <- makeTx2GeneFromGFF(file)
#' print(x)
NULL



## FIXME RETHINK THIS, CALLING AN INTERNAL GENERATOR WITHOUT THE
## BROADCLASS AND SYNONYMS OPTIONS.

#' @rdname makeTx2Gene
#' @export
## Updated 2020-10-06.
makeTx2GeneFromEnsembl <-
    function() {
        gr <- do.call(
            what = makeGRangesFromEnsembl,
            args = matchArgsToDoCall(
                args = list(
                    level = "transcripts",
                    broadClass = FALSE,
                    synonyms = FALSE
                )
            )
        )
        Tx2Gene(gr)
    }

f <- formals(makeGRangesFromEnsembl)
f <- f[setdiff(names(f), "level")]
formals(makeTx2GeneFromEnsembl) <- f



## FIXME RETHINK THIS, CALLING AN INTERNAL GENERATOR WITHOUT THE
## BROADCLASS AND SYNONYMS OPTIONS.

#' @rdname makeTx2Gene
#' @export
## Updated 2020-10-06.
makeTx2GeneFromEnsDb <- function(object) {
    gr <- do.call(
        what = makeGRangesFromEnsDb,
        args = matchArgsToDoCall(
            args = list(
                level = "transcripts",
                broadClass = FALSE,
                synonyms = FALSE
            )
        )
    )
    Tx2Gene(gr)
}

f <- formals(makeGRangesFromEnsDb)
f <- f[setdiff(names(f), "level")]
formals(makeTx2GeneFromEnsDb) <- f



## FIXME RETHINK THIS, CALLING AN INTERNAL GENERATOR WITHOUT THE
## BROADCLASS AND SYNONYMS OPTIONS.

#' @rdname makeTx2Gene
#' @export
## Updated 2020-10-06.
makeTx2GeneFromGFF <- function(file) {
    gr <- do.call(
        what = makeGRangesFromGFF,
        args = matchArgsToDoCall(args = list(
            level = "transcripts",
            broadClass = FALSE,
            synonyms = FALSE
        ))
    )
    Tx2Gene(gr)
}

f <- formals(makeGRangesFromGFF)
f <- f[setdiff(names(f), c("level", ".checkAgainstTxDb"))]
formals(makeTx2GeneFromGFF) <- f



#' @rdname makeTx2Gene
#' @export
makeTx2GeneFromGTF <- makeTx2GeneFromGFF
