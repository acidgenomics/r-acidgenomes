## FIXME ENSURE IGNOREVERSIONS WORKS HERE AS EXPECTED.



#' Make a Tx2Gene object
#'
#' @section GFF/GTF file:
#'
#' Remote URLs and compressed files are supported.
#'
#' @name makeTx2Gene
#' @note Updated 2021-01-14.
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



#' @rdname makeTx2Gene
#' @export
## Updated 2021-01-14.
makeTx2GeneFromEnsembl <-
    function() {
        gr <- do.call(
            what = .makeGRangesFromEnsembl,
            args = matchArgsToDoCall(
                args = list(
                    "broadClass" = FALSE,
                    "level" = "transcripts",
                    "synonyms" = FALSE
                )
            )
        )
        Tx2Gene(gr)
    }

f <- formals(.makeGRangesFromEnsembl)
f <- f[setdiff(names(f), c("broadClass", "level", "synonyms"))]
formals(makeTx2GeneFromEnsembl) <- f



#' @rdname makeTx2Gene
#' @export
## Updated 2021-01-14.
makeTx2GeneFromEnsDb <- function(object) {
    gr <- do.call(
        what = .makeGRangesFromEnsDb,
        args = matchArgsToDoCall(
            args = list(
                "broadClass" = FALSE,
                "level" = "transcripts",
                "synonyms" = FALSE
            )
        )
    )
    Tx2Gene(gr)
}

f <- formals(.makeGRangesFromEnsDb)
f <- f[setdiff(names(f), c("broadClass", "level", "synonyms"))]
formals(makeTx2GeneFromEnsDb) <- f



#' @rdname makeTx2Gene
#' @export
## Updated 2021-01-14.
makeTx2GeneFromGFF <- function(file) {
    gr <- do.call(
        what = .makeGRangesFromGFF,
        args = matchArgsToDoCall(
            args = list(
                "broadClass" = FALSE,
                "level" = "transcripts",
                "synonyms" = FALSE
            )
        )
    )
    Tx2Gene(gr)
}

f <- formals(.makeGRangesFromGFF)
f <- f[setdiff(names(f), c("broadClass", "level", "synonyms"))]
formals(makeTx2GeneFromGFF) <- f
