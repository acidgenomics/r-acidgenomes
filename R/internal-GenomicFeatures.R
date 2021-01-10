#' Make TxDb from GRanges
#'
#' @details
#' This step can be noisy and generate expected warnings:
#'
#' ```
#' some exons are linked to transcripts not found in the file
#' The following orphan exon were dropped
#' The following orphan CDS were dropped
#' ```
#'
#' @seealso
#' - https://stackoverflow.com/questions/38603668
#' - https://stackoverflow.com/questions/16517795
#'
#' @noRd
.makeTxDbFromGFF <- function(file) {
    alert(sprintf(
        "Making {.var %s} from {.file %s} with {.pkg %s}::{.fun %s}.",
        "TxDb", file,
        "GenomicFeatures", "makeTxDbFromGFF"
    ))
    requireNamespaces("GenomicFeatures")
    assert(isAFile(file))
    txdb <- withCallingHandlers(expr = {
        tryCatch(
            expr = {
                GenomicFeatures::makeTxDbFromGFF(file)
            },
            error = function(e) {
                stop(paste(
                    "Failed to make TxDb using ",
                    "'GenomicFeatures::makeTxDbFromGFF()'",
                    sep = "\n"
                ))
            }
        )
    }, warning = function(w) {
        if (grepl(pattern = "stop_codon", x = conditionMessage(w))) {
            invokeRestart("muffleWarning")
        }
    })
    assert(is(txdb, "TxDb"))
    validObject(txdb)
    txdb
}



## nocov end
