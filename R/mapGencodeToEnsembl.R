#' Map GENCODE release to Ensembl
#'
#' @export
#' @note Updated 2023-01-30.
#'
#' @param release `integer(1)`.
#' Human (e.g. `42`) or mouse (e.g. `"M21"`) GENCODE release.
#'
#' @return `integer(1)`.
#' Ensembl release.
#'
#' @seealso
#' - https://www.gencodegenes.org/human/releases.html
#' - https://www.gencodegenes.org/mouse/releases.html
#' - https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/_README.TXT
#' - https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/_README.TXT
#'
#' @examples
#' ## Homo sapiens.
#' mapGencodeToEnsembl(42L)
#' ## Mus musculus.
#' mapGencodeToEnsembl("M31")
mapGencodeToEnsembl <- function(release) {
    assert(isScalar(release))
    data <- import(
        con = system.file(
            "extdata", "gencode-to-ensembl.rds",
            package = .pkgName
        ),
        quiet = TRUE
    )
    assert(isSubset(c("gencode", "ensembl"), colnames(data)))
    idx <- match(x = as.character(release), table = data[["gencode"]])
    assert(
        isInt(idx),
        msg = sprintf(
            "Failed to match GENCODE release: {.var %s}.",
            release
        )
    )
    out <- data[["ensembl"]][idx]
    assert(isIntegerish(out))
    out <- as.integer(out)
    out
}
