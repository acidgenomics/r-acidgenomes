#' Get the AnnotationHub identifier for desired EnsDb
#'
#' @note Updated 2023-12-04.
#' @noRd
#'
#' @examples
#' .getEnsDbAnnotationHubId("Homo sapiens")
#' .getEnsDbAnnotationHubId("Canis lupus familiaris")
.getEnsDbAnnotationHubId <-
    function(organism, genomeBuild = NULL, release = NULL, ah = NULL) {
        assert(
            isOrganism(organism),
            isString(genomeBuild, nullOk = TRUE),
            isInt(release, nullOk = TRUE),
            is(ah, "AnnotationHub") || is.null(ah)
        )
        ## Coerce integerish value (e.g. "90") to integer (e.g. "90L").
        if (isInt(release)) {
            release <- as.integer(release)
        }
        ## Get AnnotationHub.
        if (is.null(ah)) {
            ah <- .annotationHub()
        }
        ## Matching EnsDb objects from ensembldb by default.
        preparerClass <- "AHEnsDbs"
        rdataClass <- "EnsDb"
        alert(sprintf(
            "Getting {.cls %s} from {.pkg %s} %s (%s).",
            rdataClass,
            "AnnotationHub",
            packageVersion("AnnotationHub"),
            AnnotationHub::snapshotDate(ah)
        ))
        ## Query AnnotationHub.
        ahs <- AnnotationHub::query(
            x = ah,
            pattern = c(
                "Ensembl",
                organism,
                genomeBuild,
                preparerClass,
                rdataClass,
                release
            ),
            ignore.case = TRUE
        )
        assert(is(ahs, "AnnotationHub"))
        ## Get the AnnotationHub from the metadata columns.
        mcols <- mcols(ahs, use.names = TRUE)
        assert(isSubset(
            x = c(
                "dataprovider",
                "genome",
                "species",
                "tags",
                "title"
            ),
            y = colnames(mcols)
        ))
        keep <- grepl(
            pattern = paste0("^", organism),
            x = mcols[["species"]]
        )
        mcols <- mcols[keep, , drop = FALSE]
        ## Ensure genome build matches exactly.
        if (!is.null(genomeBuild)) {
            keep <- mcols[["genome"]] == genomeBuild
            mcols <- mcols[keep, , drop = FALSE]
        }
        ## Ensure release version matches exactly.
        if (!is.null(release)) {
            keep <- grepl(
                pattern = paste0("\\b", release, "\\b"),
                x = mcols[["title"]]
            )
            mcols <- mcols[keep, , drop = FALSE]
            keep <- vapply(
                X = mcols[["tags"]],
                FUN = function(tags) {
                    release %in% tags
                },
                FUN.VALUE = logical(1L),
                USE.NAMES = FALSE
            )
            mcols <- mcols[keep, , drop = FALSE]
        }
        assert(
            hasRows(mcols),
            all(mcols[["dataprovider"]] == "Ensembl"),
            all(mcols[["preparerclass"]] == preparerClass),
            all(mcols[["rdataclass"]] == rdataClass),
            all(mcols[["sourcetype"]] == "ensembl"),
            msg = "No entry matched on AnnotationHub."
        )
        ## Sort the entries by Ensembl release as integer instead of
        ## AnnotationHub identifier. Updates can otherwise mess up the expected
        ## order, for example:
        ## > AH73881 | Ensembl 97 EnsDb for Homo sapiens
        ## > AH73986 | Ensembl 79 EnsDb for Homo sapiens
        ## > AH79689 | Ensembl 100 EnsDb for Homo sapiens
        match <- strMatch(
            x = mcols[["title"]],
            pattern = "^Ensembl ([0-9]+) EnsDb.+$"
        )
        idx <- order(as.integer(match[, 2L]))
        mcols <- mcols[idx, , drop = FALSE]
        ## Error if filtering was unsuccessful.
        if (!hasRows(mcols)) {
            abort(sprintf(
                fmt = paste(
                    "No entry matched on AnnotationHub {.val %s}.",
                    "  - {.arg %s}: {.val %s}",
                    "  - {.arg %s}: {.val %s}",
                    "  - {.arg %s}: {.val %s}",
                    sep = "\n"
                ),
                as.character(packageVersion("AnnotationHub")),
                "Organism",
                as.character(organism),
                "Genome build",
                as.character(genomeBuild),
                "Ensembl release",
                as.character(release)
            ))
        }
        ## Select the most recent database (sorted by title, not identifier!).
        mcols <- tail(mcols, n = 1L)
        id <- rownames(mcols)
        assert(
            isString(id),
            isMatchingRegex(x = id, pattern = "^AH[[:digit:]]+$")
        )
        alertInfo(sprintf("{.val %s}: %s.", id, mcols[["title"]]))
        id
    }


#' Get EnsDb from AnnotationHub identifier
#'
#' @note Updated 2023-04-26.
#' @noRd
#'
#' @details
#' This step will also output `txProgressBar` on a fresh install. Using
#' `capture.output` here again to suppress console output. Additionally, it
#' attaches ensembldb and other Bioconductor dependency packages, which will
#' mask some tidyverse functions (e.g. `select`).
#'
#' @examples
#' edb <- .getEnsDbFromAnnotationHub("AH64923")
#' print(edb)
.getEnsDbFromAnnotationHub <- function(id, ah = NULL) {
    assert(
        isString(id),
        is(ah, "AnnotationHub") || is.null(ah)
    )
    if (is.null(ah)) {
        ah <- .annotationHub()
    }
    assert(is(ah, "AnnotationHub"))
    quietly({
        edb <- ah[[id]]
    })
    assert(is(edb, "EnsDb"))
    edb
}


#' Get EnsDb from Package
#'
#' @note Updated 2023-04-26.
#' @noRd
#'
#' @examples .getEnsDbFromPackage("EnsDb.Hsapiens.v75")
.getEnsDbFromPackage <- function(package) {
    alert(sprintf("Getting {.cls %s} from {.pkg %s}.", "EnsDb", package))
    assert(
        isString(package),
        requireNamespaces(package)
    )
    edb <- get(
        x = package,
        envir = asNamespace(package),
        inherits = FALSE
    )
    assert(is(edb, "EnsDb"))
    edb
}


#' Get metadata inside EnsDb object
#'
#' @note Updated 2023-04-26.
#' @noRd
.getEnsDbMetadata <- function(object, level = NULL) {
    assert(
        requireNamespaces("ensembldb"),
        is(object, "EnsDb"),
        isString(level, nullOk = TRUE)
    )
    metadata <- metadata(object)
    assert(is.data.frame(metadata))
    genomeBuild <- metadata[
        match(x = "genome_build", table = metadata[["name"]]),
        "value",
        drop = TRUE
    ]
    assert(isString(genomeBuild))
    organism <- organism(object)
    organism <- switch(
        EXPR = organism,
        "Saccharomyces cerevisiae S288c" = {
            ## AH109740 Ensembl 109.
            "Saccharomyces cerevisiae"
        },
        organism
    )
    quietly({
        release <- ensembldb::ensemblVersion(object)
    })
    release <- as.integer(release)
    list <- list(
        "ensembldb" = metadata,
        "genomeBuild" = genomeBuild,
        "organism" = organism,
        "provider" = "Ensembl",
        "release" = release
    )
    if (!is.null(level)) {
        list[["level"]] <- level
    }
    ## AnnotationHub identifier should be stashed in attributes, when possible.
    if (isString(attr(object, "annotationHubId"))) {
        list[["annotationHubId"]] <- attr(object, "annotationHubId")
    }
    items <- c(
        "Organism" = list[["organism"]],
        "Genome build" = list[["genomeBuild"]],
        "Release" = list[["release"]]
    )
    if (isString(list[["level"]])) {
        items <- c(items, "Level" = list[["level"]])
    }
    dl(items)
    list
}
