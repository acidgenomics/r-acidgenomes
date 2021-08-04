#' Get EnsDb from Bioconductor
#'
#' @export
#' @note Updated 2021-08-04.
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
#' edb <- getEnsDb(organism = "Homo sapiens", release = 100L)
#' print(edb)
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
    organism <- gsub(pattern = "_", replacement = " ", x = makeNames(organism))
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



## FIXME Should we not sanitize the genome Build to simple name here?

#' Get the AnnotationHub identifier for desired EnsDb
#'
#' @note Updated 2021-08-03.
#' @noRd
#'
#' @examples
#' .getAnnotationHubID("Homo sapiens")
.getEnsDbAnnotationHubID <- function(
    organism,
    genomeBuild = NULL,
    release = NULL,
    ah = NULL
) {
    requireNamespaces("AnnotationHub")
    assert(
        isString(organism),
        isString(genomeBuild, nullOK = TRUE),
        isInt(release, nullOK = TRUE),
        is(ah, "AnnotationHub") || is.null(ah)
    )
    ## Standardize organism name, if necessary.
    organism <- gsub(pattern = "_", replacement = " ", x = makeNames(organism))
    ## ensembldb always uses two words for organisms, instead of matching the
    ## Ensembl name exactly. This can mismatch with some organisms. For example,
    ## the dog genome is named "Canis lupus familiaris" on Ensembl but matches
    ## against "Canis familiaris" only with ensembldb. Check for this rare edge
    ## case and inform the user.
    pattern <- "^([a-z]+)\\s[a-z]+\\s([a-z]+)$"
    if (isTRUE(grepl(pattern = pattern, x = organism, ignore.case = TRUE))) {
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
    ## Coerce integerish (e.g. 90) to integer (e.g. 90L).
    if (isInt(release)) {
        release <- as.integer(release)
    }
    ## Error on request of unsupported legacy Ensembl release.
    ## Don't hardcode against Ensembl 87 cutoff, in case older releases are
    ## added back in a future ensembldb/AnnotationHub update.
    ## > if (
    ## >     is.integer(release) &&
    ## >     release < 87L
    ## > ) {
    ## >     stop("ensembldb currently only supports Ensembl releases >= 87.")
    ## > }
    ## Get AnnotationHub.
    if (is.null(ah)) {
        ah <- .annotationHub()
    }
    ## Matching EnsDb objects from ensembldb by default.
    preparerclass <- "AHEnsDbs"
    rdataclass <- "EnsDb"
    alert(sprintf(
        "Getting {.var %s} from {.pkg %s} %s (%s).",
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
    ## FIXME This step is failing for "Drosophila melanogaster" with "BDGP6".
    ## FIXME This needs to match the genomeBuild more precisely.
    ## "BDGP6" needs to match exactly and not allow propagation of "BDGP6.32",
    ## for example.
    assert(
        all(mcols[["dataprovider"]] == "Ensembl"),
        all(mcols[["preparerclass"]] == preparerclass),
        all(mcols[["rdataclass"]] == rdataclass),
        all(mcols[["sourcetype"]] == "ensembl"),
        all(tolower(mcols[["genome"]]) == tolower(genomeBuild)),
        all(tolower(mcols[["species"]]) == tolower(organism)),
        msg = "Invalid metadata returned from AnnotationHub query."
    )
    ## Sort the entries by Ensembl release as integer instead of AH identifier.
    ## Updates can otherwise mess up the expected order, for example:
    ## > AH73881 | Ensembl 97 EnsDb for Homo sapiens
    ## > AH73986 | Ensembl 79 EnsDb for Homo sapiens
    ## > AH79689 | Ensembl 100 EnsDb for Homo sapiens
    match <- str_match(
        string = mcols[["title"]],
        pattern = "^Ensembl ([0-9]+) EnsDb.+$"
    )
    idx <- order(as.integer(match[, 2L]))
    mcols <- mcols[idx, , drop = FALSE]
    ## Abort if there's no match and working offline.
    if (!isTRUE(hasInternet()) && !hasRows(mcols)) {
        ## nocov start
        stop("AnnotationHub requires an Internet connection.")
        ## nocov end
    }
    ## Ensure genome build matches, if specified.
    if (!is.null(genomeBuild)) {
        assert(isSubset("genome", colnames(mcols)))
        keep <- which(mcols[["genome"]] %in% genomeBuild)
        mcols <- mcols[keep, , drop = FALSE]
    }
    ## Ensure Ensembl release matches, or pick the latest one.
    if (!is.null(release)) {
        assert(isSubset("title", colnames(mcols)))
        keep <- which(grepl(paste("Ensembl", release), mcols[["title"]]))
        mcols <- mcols[keep, , drop = FALSE]
        assert(hasLength(nrow(mcols), n = 1L))
    }
    ## Error if filtering was unsuccessful.
    if (!hasRows(mcols)) {
        stop(sprintf(
            fmt = paste(
                "No entry matched on AnnotationHub %s.",
                "  - %s: %s",
                "  - %s: %s",
                "  - %s: %s",
                sep = "\n"
            ),
            packageVersion("AnnotationHub"),
            "Organism", deparse(organism),
            "Genome build", deparse(genomeBuild),
            "Ensembl release", deparse(release)
        ))
    }
    ## Select the most recent database (sorted by title, not identifier!).
    mcols <- tail(mcols, n = 1L)
    id <- rownames(mcols)
    assert(
        isString(id),
        unname(isMatchingRegex(x = id, pattern = "^AH[[:digit:]]+$"))
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
#' @note Updated 2021-01-18.
#' @noRd
#'
#' @examples .getEnsDbFromPackage("EnsDb.Hsapiens.v75")
.getEnsDbFromPackage <- function(package) {
    alert(sprintf("Getting {.var %s} from {.pkg %s}.", "EnsDb", package))
    assert(isString(package))
    require(package, character.only = TRUE)
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
