#' Get EnsDb from Bioconductor AnnotationHub
#'
#' @export
#' @note Updated 2020-09-25.
#'
#' @inheritParams AcidRoxygen::params
#'
#' @details
#' Remaps UCSC genome build to Ensembl automatically, if necessary.
#' Provides legacy support for GRCh37 (hg19).
#'
#' @return `EnsDb`.
#'
#' @examples
#' getEnsDb(organism = "Homo sapiens", release = 100L)
getEnsDb <- function(
    organism,
    genomeBuild = NULL,
    release = NULL
) {
    assert(
        isString(organism),
        isString(genomeBuild, nullOK = TRUE),
        isInt(release, nullOK = TRUE)
    )
    if (isString(genomeBuild)) {
        remap <- tryCatch(
            expr = convertUCSCBuildToEnsembl(genomeBuild),
            error = function(e) NULL
        )
        if (hasLength(remap)) {
            ucsc <- names(remap)
            ensembl <- unname(remap)
            cli_alert_warning(sprintf(
                fmt = paste(
                    "Remapping genome build from UCSC ({.val %s}) to",
                    "Ensembl ({.val %s})."
                ),
                ucsc, ensembl
            ))
            genomeBuild <- ensembl
            rm(remap, ucsc, ensembl)
        }
    }
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
        id <- .getAnnotationHubID(
            organism = organism,
            genomeBuild = genomeBuild,
            release = release
        )
        edb <- .getEnsDbFromAnnotationHub(id = id)
    }
    attr(edb, "id") <- id
    edb
}



#' Get EnsDb from AnnotationHub
#'
#' @note Updated 2020-09-24.
#' @noRd
#'
#' @details
#' This step will also output `txProgressBar` on a fresh install. Using
#' `capture.output` here again to suppress console output. Additionally, it
#' attaches ensembldb and other Bioconductor dependency packages, which will
#' mask some tidyverse functions (e.g. `select`).
#'
#' @examples
#' .getEnsDbFromAnnotationHub("AH64923")
.getEnsDbFromAnnotationHub <- function(id, ah = NULL) {
    assert(
        isString(id),
        is(ah, "AnnotationHub") || is.null(ah)
    )
    userAttached <- .packages()
    if (is.null(ah)) {
        ah <- .annotationHub()
    }
    assert(is(ah, "AnnotationHub"))
    invisible(capture.output({
        edb <- suppressMessages(ah[[id]])
    }))
    assert(is(edb, "EnsDb"))
    .forceDetach(keep = userAttached)
    edb
}



#' Get EnsDb from Package
#'
#' @note Updated 2020-01-20.
#' @noRd
#'
#' @examples .getEnsDbFromPackage("EnsDb.Hsapiens.v75")
.getEnsDbFromPackage <- function(package) {
    cli_alert(sprintf("Getting {.var EnsDb} from {.pkg %s}.", package))
    userAttached <- .packages()
    assert(isString(package))
    require(package, character.only = TRUE)
    edb <- get(
        x = package,
        envir = asNamespace(package),
        inherits = FALSE
    )
    assert(is(edb, "EnsDb"))
    .forceDetach(keep = userAttached)
    edb
}



#' Get metadata inside EnsDb object
#'
#' @note Updated 2020-09-25.
#' @noRd
.getEnsDbMetadata <- function(object, level = NULL) {
    assert(
        is(object, "EnsDb"),
        isString(level, nullOK = TRUE)
    )
    metadata <- metadata(object)
    assert(is.data.frame(metadata))
    genomeBuild <- metadata[
        match(x = "genome_build", table = metadata[["name"]]),
        "value",
        drop = TRUE
    ]
    assert(isString(genomeBuild))
    list <- list(
        "organism" = organism(object),
        "genomeBuild" = genomeBuild,
        "ensemblRelease" = as.integer(ensemblVersion(object)),
        "ensembldb" = metadata
    )
    if (!is.null(level)) {
        level[["level"]] <- level
    }
    ## AnnotationHub ID should be stashed in attributes, when possible.
    if (isString(attr(object, "id"))) {
        list[["id"]] <- attr(object, "id")
    }
    cli_div(theme = list(body = list("margin-left" = 4L)))
    items <- c(
        "Organism" = list[["organism"]],
        "Genome build" = list[["genomeBuild"]],
        "Ensembl release" = list[["ensemblRelease"]]
    )
    if (isString(list[["level"]])) {
        items <- c(items, "Level" = list[["level"]])
    }
    cli_dl(items = items)
    cli_end()
    list
}
