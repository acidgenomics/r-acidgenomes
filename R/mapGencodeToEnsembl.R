#' Map GENCODE release to Ensembl
#'
#' @export
#' @note Updated 2022-01-17.
#'
#' @param release `integer(1)`.
#'   Human GENCODE release (e.g. `39`).
#'
#' @return `integer(1)`.
#'   Ensembl release.
#'
#' @seealso
#' - https://www.gencodegenes.org/
#' - https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/_README.TXT
#'
#' @examples
#' mapGencodeToEnsembl(39L)
mapGencodeToEnsembl <- function(release) {
    assert(isInt(release))
    data <- import(
        file = system.file(
            "extdata", "gencode-to-ensembl.rds",
            package = .pkgName
        ),
        quiet = TRUE
    )
    assert(isSubset(c("gencode", "ensembl"), colnames(data)))
    idx <- match(x = release, table = data[["gencode"]])
    assert(
        isInt(idx),
        msg = "Failed to match GENCODE release."
    )
    out <- data[["ensembl"]][idx]
    assert(isInt(out))
    out <- as.integer(out)
    out
}
