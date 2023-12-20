#' Get EnsDb from Bioconductor
#'
#' @export
#' @note Updated 2023-12-20.
#'
#' @inheritParams AcidRoxygen::params
#'
#' @details
#' Primarily fetches `EnsDb` objects from AnnotationHub.
#' For legacy GRCh37 genome build, loads `EnsDb.Hsapiens.v75` data package.
#'
#' @section UCSC genome build remapping:
#'
#' Remaps UCSC genome build to Ensembl automatically, if necessary.
#' For example, `"hg38"` remaps to `"GRCh38"` and "hg19" to `"GRCh37"`.
#'
#' @return `EnsDb`.
#'
#' @examples
#' ## Current genome build and release.
#' object <- getEnsDb(organism = "Homo sapiens")
#' print(object)
#'
#' ## Legacy GRCh37 genome build.
#' if (requireNamespace("EnsDb.Hsapiens.v75", quietly = TRUE)) {
#'     object <- getEnsDb(organism = "Homo sapiens", genomeBuild = "GRCh37")
#'     print(object)
#' }
getEnsDb <-
    function(organism,
             genomeBuild = NULL,
             release = NULL) {
        assert(
            isString(organism),
            isString(genomeBuild, nullOk = TRUE),
            isInt(release, nullOk = TRUE)
        )
        organism <- gsub(
            pattern = "_",
            replacement = " ",
            x = makeNames(organism)
        )
        if (
            identical(tolower(organism), "homo sapiens") &&
            (
                identical(tolower(as.character(genomeBuild)), "grch37") ||
                identical(release, 75L)
            )
        ) {
            id <- "EnsDb.Hsapiens.v75"
            edb <- .getEnsDbFromPackage(package = id)
        } else {
            id <- .getEnsDbAnnotationHubId(
                organism = organism,
                genomeBuild = genomeBuild,
                release = release
            )
            edb <- .getEnsDbFromAnnotationHub(id = id)
        }
        attr(edb, "annotationHubId") <- id
        edb
    }
