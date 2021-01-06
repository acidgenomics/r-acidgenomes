#' Connect to AnnotationHub
#'
#' @note Updated 2020-09-24.
#' @noRd
#'
#' @details
#' On a fresh install this will print a txProgressBar to the console. We're
#' using [utils::capture.output()] here to suppress the console output, since
#' it's not very informative and can cluster R Markdown reports.
.annotationHub <- function() {
    userAttached <- .packages()
    invisible(capture.output({
        ah <- suppressMessages(AnnotationHub())
    }))
    assert(is(ah, "AnnotationHub"))
    forceDetach(keep = userAttached)
    ah
}



#' Get the AnnotationHub ID for desired EnsDb
#'
#' @note Updated 2020-10-12.
#' @noRd
#'
#' @examples
#' .getAnnotationHubID("Homo sapiens")
.getAnnotationHubID <- function(
    organism,
    genomeBuild = NULL,
    release = NULL,
    ah = NULL
) {
    assert(
        isOrganism(organism),
        isString(genomeBuild, nullOK = TRUE),
        isInt(release, nullOK = TRUE),
        is(ah, "AnnotationHub") || is.null(ah)
    )
    userAttached <- .packages()
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
        cli_alert(sprintf(
            "Matching {.val %s} using {.val %s}.",
            fullOrganism, organism
        ))
    }
    ## Coerce integerish (e.g. 90) to integer (e.g. 90L).
    if (isInt(release)) {
        release <- as.integer(release)
    }
    ## Error on request of unsupported legacy Ensembl release.
    if (
        is.integer(release) &&
        release < 87L
    ) {
        stop("ensembldb currently only supports Ensembl releases >= 87.")
    }
    ## Get AnnotationHub.
    if (is.null(ah)) {
        ah <- .annotationHub()
    }
    ## Matching EnsDb objects from ensembldb by default.
    preparerclass <- "AHEnsDbs"
    rdataclass <- "EnsDb"
    cli_alert(sprintf(
        "Getting {.var %s} from {.pkg AnnotationHub} %s (%s).",
        rdataclass,
        packageVersion("AnnotationHub"),
        snapshotDate(ah)
    ))
    ## Query AnnotationHub.
    ahs <- query(
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
    assert(
        all(mcols[["dataprovider"]] == "Ensembl"),
        all(mcols[["preparerclass"]] == preparerclass),
        all(mcols[["rdataclass"]] == rdataclass),
        all(mcols[["sourcetype"]] == "ensembl"),
        all(tolower(mcols[["genome"]]) == tolower(genomeBuild)),
        all(tolower(mcols[["species"]]) == tolower(organism))
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
                "No ID matched on AnnotationHub %s.",
                "  - Organism: %s",
                "  - Genome build: %s",
                "  - Ensembl release: %s",
                sep = "\n"
            ),
            packageVersion("AnnotationHub"),
            deparse(organism),
            deparse(genomeBuild),
            deparse(release)
        ))
    }
    ## Select the most recent database (sorted by title, not ID!).
    mcols <- tail(mcols, n = 1L)
    id <- rownames(mcols)
    assert(
        isString(id),
        unname(isMatchingRegex(x = id, pattern = "^AH[[:digit:]]+$"))
    )
    cli_alert_info(sprintf("{.val %s}: %s.", id, mcols[["title"]]))
    forceDetach(keep = userAttached)
    id
}



#' Get Ensembl/Entrez mappings from NCBI OrgDb via AnnotationHub
#'
#' @note Updated 2020-10-01.
#' @noRd
.getEnsembl2EntrezFromOrgDb <- function(
    keys,
    keytype,
    columns,
    organism,
    strict = TRUE
) {
    assert(
        isCharacter(keys),
        hasNoDuplicates(keys),
        isString(keytype),
        isCharacter(columns),
        isOrganism(organism),
        isFlag(strict)
    )
    cli_alert(sprintf(
        "Matching identifiers using NCBI {.var %s} via {.pkg %s} %s.",
        "OrgDb",
        "AnnotationHub",
        packageVersion("AnnotationHub")
    ))
    userAttached <- .packages()
    ah <- .annotationHub()
    ahs <- query(ah, pattern = c(organism, "NCBI", "OrgDb"))
    id <- tail(names(ahs), n = 1L)
    suppressMessages({
        orgdb <- ah[[id]]
    })
    assert(is(orgdb, "OrgDb"))
    cli_alert_info(sprintf(
        "{.val %s} (%s): %s.",
        id,
        mcols(ahs)[id, "title"],
        mcols(ahs)[id, "description"]
    ))
    suppressMessages({
        df <- select(
            x = orgdb,
            keys = keys,
            keytype = keytype,
            columns = columns
        )
    })
    df <- as(df, "DataFrame")
    if (isTRUE(strict)) {
        df <- df[complete.cases(df), ]
        if (!areSetEqual(keys, unique(df[[keytype]]))) {
            setdiff <- setdiff(keys, unique(df[[keytype]]))
            stop(sprintf(
                "Match failure: %s.",
                toString(setdiff, width = 200L)
            ))
        }
    }
    colnames(df)[colnames(df) == "ENSEMBL"] <- "ensembl"
    colnames(df)[colnames(df) == "ENTREZID"] <- "entrez"
    df[["entrez"]] <- as.integer(df[["entrez"]])
    forceDetach(keep = userAttached)
    df
}
