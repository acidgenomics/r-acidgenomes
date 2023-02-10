#' Map GENCODE release to Ensembl
#'
#' @export
#' @note Updated 2023-02-10.
#'
#' @param release `integer(1)` or `character(1)`.
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
    organism <- ifelse(
        test = grepl(pattern = "^M", x = release),
        yes = "Mus musculus",
        no = "Homo sapiens"
    )
    df <- gencodeReleaseHistory(organism = organism)
    assert(isSubset(c("gencodeRelease", "ensemblRelease"), colnames(df)))
    idx <- match(x = as.character(release), table = df[["gencodeRelease"]])
    assert(
        isInt(idx),
        msg = sprintf(
            "Failed to match GENCODE release: {.var %s}.",
            release
        )
    )
    out <- as.numeric(df[["ensemblRelease"]][idx])
    assert(isIntegerish(out))
    out <- as.integer(out)
    out
}
