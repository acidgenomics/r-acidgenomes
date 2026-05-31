#' Map Ensembl release to GENCODE
#'
#' @export
#' @note Updated 2026-05-31.
#'
#' @param release `integer(1)`.
#' Ensembl release (e.g. `110`).
#'
#' @param organism `character(1)`.
#' Scientific name of the organism.
#' Currently only `"Homo sapiens"` and `"Mus musculus"` are supported.
#'
#' @return `character(1)`.
#' Corresponding GENCODE release identifier (e.g. `"44"` for human or
#' `"M33"` for mouse).
#'
#' @seealso
#' - https://www.gencodegenes.org/human/releases.html
#' - https://www.gencodegenes.org/mouse/releases.html
#' - [mapGencodeToEnsembl()].
#'
#' @examples
#' ## Homo sapiens.
#' mapEnsemblToGencode(110L, organism = "Homo sapiens")
#' ## Mus musculus.
#' mapEnsemblToGencode(110L, organism = "Mus musculus")
mapEnsemblToGencode <- function(release, organism = "Homo sapiens") {
    assert(
        isInt(release),
        isOrganism(organism)
    )
    organism <- match.arg(
        arg = organism,
        choices = c("Homo sapiens", "Mus musculus")
    )
    df <- gencodeReleaseHistory(organism = organism)
    assert(isSubset(c("gencodeRelease", "ensemblRelease"), colnames(df)))
    idx <- match(
        x = as.character(as.integer(release)),
        table = as.character(as.numeric(df[["ensemblRelease"]]))
    )
    if (is.na(idx)) {
        abort(sprintf(
            "Failed to match Ensembl release {.val %s} for {.emph %s}.",
            release,
            organism
        ))
    }
    as.character(df[["gencodeRelease"]][idx])
}
