#' Get EnsDb from Bioconductor
#'
#' @note Updated 2022-01-12.
#' @noRd
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
#' edb <- .getEnsDb(organism = "Homo sapiens", release = 100L)
#' print(edb)
.getEnsDb <-
    function(organism,
             genomeBuild = NULL,
             release = NULL) {
        assert(
            isString(organism),
            isString(genomeBuild, nullOK = TRUE),
            isInt(release, nullOK = TRUE)
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
            id <- .getEnsDbAnnotationHubID(
                organism = organism,
                genomeBuild = genomeBuild,
                release = release
            )
            edb <- .getEnsDbFromAnnotationHub(id = id)
        }
        attr(edb, "annotationHubId") <- id
        edb
    }



#' Get the AnnotationHub identifier for desired EnsDb
#'
#' @note Updated 2021-08-04.
#' @noRd
#'
#' @examples
#' .getAnnotationHubID("Homo sapiens")
.getEnsDbAnnotationHubID <-
    function(organism,
             genomeBuild = NULL,
             release = NULL,
             ah = NULL) {
        assert(
            requireNamespaces("AnnotationHub"),
            isString(organism),
            isString(genomeBuild, nullOK = TRUE),
            isInt(release, nullOK = TRUE),
            is(ah, "AnnotationHub") || is.null(ah)
        )
        ## Standardize organism name, if necessary.
        organism <- gsub(
            pattern = "_",
            replacement = " ",
            x = makeNames(organism)
        )
        ## The ensembldb package always uses two words for organisms, instead of
        ## matching the Ensembl name exactly. This can mismatch with some
        ## organisms. For example, the dog genome is named "Canis lupus
        ## familiaris" on Ensembl but matches against "Canis familiaris" only
        ## with ensembldb. Check for this rare edge case and inform the user.
        pattern <- "^([a-z]+)\\s[a-z]+\\s([a-z]+)$"
        if (isTRUE(grepl(
            pattern = pattern,
            x = organism,
            ignore.case = TRUE
        ))) {
            fullOrganism <- organism
            organism <- sub(
                pattern = pattern,
                replacement = "\\1 \\2",
                x = fullOrganism,
                ignore.case = TRUE
            )
            alert(sprintf(
                "Matching {.val %s} using {.val %s}.",
                fullOrganism, organism
            ))
        }
        assert(isOrganism(organism))
        ## Coerce integerish value (e.g. "90") to integer (e.g. "90L").
        if (isInt(release)) {
            release <- as.integer(release)
        }
        ## Get AnnotationHub.
        if (is.null(ah)) {
            ah <- .annotationHub()
        }
        ## Matching EnsDb objects from ensembldb by default.
        preparerclass <- "AHEnsDbs"
        rdataclass <- "EnsDb"
        alert(sprintf(
            "Getting {.cls %s} from {.pkg %s} %s (%s).",
            rdataclass,
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
                preparerclass,
                rdataclass,
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
        ## Ensure organism (species) matches exactly.
        keep <- mcols[["species"]] == organism
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
            all(mcols[["preparerclass"]] == preparerclass),
            all(mcols[["rdataclass"]] == rdataclass),
            all(mcols[["sourcetype"]] == "ensembl"),
            msg = "No entry matched on AnnotationHub."
        )
        ## Sort the entries by Ensembl release as integer instead of
        ## AnnotationHub identifier. Updates can otherwise mess up the expected
        ## order, for example:
        ## > AH73881 | Ensembl 97 EnsDb for Homo sapiens
        ## > AH73986 | Ensembl 79 EnsDb for Homo sapiens
        ## > AH79689 | Ensembl 100 EnsDb for Homo sapiens
        match <- stri_match_first_regex(
            str = mcols[["title"]],
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
                "Organism", as.character(organism),
                "Genome build", as.character(genomeBuild),
                "Ensembl release", as.character(release)
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
#' @note Updated 2021-01-14.
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
    invisible(capture.output({
        edb <- suppressMessages(ah[[id]])
    }))
    assert(is(edb, "EnsDb"))
    edb
}



#' Get EnsDb from Package
#'
#' @note Updated 2021-09-13.
#' @noRd
#'
#' @examples .getEnsDbFromPackage("EnsDb.Hsapiens.v75")
.getEnsDbFromPackage <- function(package) {
    alert(sprintf("Getting {.cls %s} from {.pkg %s}.", "EnsDb", package))
    assert(isString(package))
    requireNamespaces(package)
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
#' @note Updated 2021-04-27.
#' @noRd
.getEnsDbMetadata <- function(object, level = NULL) {
    requireNamespaces("ensembldb")
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
        "ensembldb" = metadata,
        "genomeBuild" = genomeBuild,
        "organism" = organism(object),
        "provider" = "Ensembl",
        "release" = as.integer(ensembldb::ensemblVersion(object))
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
