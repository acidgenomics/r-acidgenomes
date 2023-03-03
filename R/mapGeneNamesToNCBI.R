## FIXME Add support for this.
## FIXME Need to add support for



#' Map gene names to NCBI
#'
#' @export
#' @note Updated 2023-03-03.
#'
#' @inheritParams AcidRoxygen::params
#'
#' @param genes `character`.
#' Gene names (e.g. `"TUT4"`).
#'
#' @param ncbi `NcbiGeneIfo` or `NULL`.
#' If `NULL`, NCBI annotations will be downloaded automatically.
#'
#' @examples
#' x <- mapGeneNamesToNCBI(
#'     genes = c("TUT4", "ZCCHC11", "TENT3A"),
#'     organism = "Homo sapiens"
#' )
mapGeneNamesToNCBI <-
    function(
        genes,
        organism,
        ncbi
    ) {
        ## FIXME Ensure that organism and NCBI match.

    }
