#' Get the directives from a GFF file
#'
#' @note Updated 2023-11-29.
#' @noRd
#'
#' @inheritParams AcidRoxygen::params
#'
#' @details
#' Matches lines beginning with `#!<key> <value>` or `##<key>: <value>`
#'
#' @section GFF3:
#'
#' Lines beginning with '##' are directives (sometimes called pragmas or
#' meta-data) and provide meta-information about the document as a whole. Blank
#' lines should be ignored by parsers and lines beginning with a single '#' are
#' used for human-readable comments and can be ignored by parsers. End-of-line
#' comments (comments preceded by # at the end of and on the same line as a
#' feature or directive line) are not allowed.
#'
#' @seealso
#' - https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
#'
#' @return `DFrame` or `NULL`.
#'
#' @examples
#' url <- pasteUrl(
#'     "ftp.ensembl.org",
#'     "pub",
#'     "release-102",
#'     "gtf",
#'     "homo_sapiens",
#'     "Homo_sapiens.GRCh38.102.gtf.gz",
#'     protocol = "https"
#' )
#' df <- .getGffDirectives(url)
#' print(df)
.getGffDirectives <- function(file, nMax = Inf) {
    assert(.isSupportedGff(file))
    file <- .cacheIt(file)
    lines <- import(
        con = file,
        format = "lines",
        comment = "",
        nMax = nMax,
        quiet = TRUE
    )
    pattern <- "^(#!|#+)([a-z-]+)(:)?\\s+(.+)$"
    lines <- grep(pattern = pattern, x = lines, value = TRUE)
    if (!hasLength(lines)) {
        return(NULL)
    }
    mat <- strMatch(x = lines, pattern = pattern)
    assert(is.matrix(mat), hasRows(mat))
    df <- as(mat, "DFrame")
    df <- df[, c(3L, 5L), drop = FALSE]
    colnames(df) <- c("key", "value")
    df <- unique(df)
    df <- df[order(df[["key"]]), , drop = FALSE]
    ## GENCODE GRCh37 GTF file incorrectly returns format as "gff3" currently.
    ## Filing issue with gencode-help@ebi.ac.uk.
    if (
        isSubset(c("format", "provider"), df[["key"]]) &&
            identical(
                x = df[["value"]][which(df[["key"]] == "format")],
                y = "gff3"
            ) &&
            identical(
                x = df[["value"]][which(df[["key"]] == "provider")],
                y = "GENCODE"
            ) &&
            isMatchingFixed(x = fileExt(file), pattern = "gtf") &&
            isMatchingFixed(
                x = df[["value"]][which(df[["key"]] == "description")],
                pattern = "GRCh37"
            )
    ) {
        df[["value"]][which(df[["key"]] == "format")] <- "gtf"
    }
    df
}



## FIXME This is incorrectly returning FlyBase for Ensembl file:
## https://ftp.ensembl.org/pub/release-110/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.46.110.gtf.gz

## FIXME Need to test detection support of:
## https://ftp.ensembl.org/pub/release-110/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.110.gtf.gz

#' Get metadata about a GFF file
#'
#' @note Updated 2024-01-02.
#' @noRd
#'
#' @inheritParams AcidRoxygen::params
#'
#' @return `list`.
#' Containing values, when possible:
#' - `directives`
#' - `format`
#' - `genomeBuild`
#' - `gffVersion`
#' - `md5`
#' - `organism`
#' - `provider`
#' - `release`
#' - `sha256`
#'
#' @examples
#' url <- pasteUrl(
#'     "ftp.ensembl.org",
#'     "pub",
#'     "release-102",
#'     "gtf",
#'     "homo_sapiens",
#'     "Homo_sapiens.GRCh38.102.gtf.gz",
#'     protocol = "https"
#' )
#' x <- .getGffMetadata(url)
#' print(x)
.getGffMetadata <- function(file) {
    assert(.isSupportedGff(file))
    file <- .cacheIt(file)
    l <- list()
    if (isAFile(file)) {
        file <- realpath(file)
    }
    l[["file"]] <- file
    l[["md5"]] <- .md5(file)
    l[["sha256"]] <- .sha256(file)
    ## Attempt to get genome build and provider from GFF directives.
    df <- .getGffDirectives(file)
    if (is(df, "DFrame")) {
        l[["directives"]] <- df
        ## These are GFF specific (not defined in GTF), but useful:
        l[["gffVersion"]] <-
            df[df[["key"]] == "gff-version", "value", drop = TRUE]
        l[["format"]] <-
            toupper(df[df[["key"]] == "format", "value", drop = TRUE])
        l[["genomeBuild"]] <- .gffGenomeBuild(df)
        l[["provider"]] <- .gffProvider(df)
    }
    if (!isString(l[["format"]])) {
        l[["format"]] <- .gffFormat(file)
    }
    assert(isSubset(l[["format"]], c("GFF3", "GTF")))
    if (isString(l[["gffVersion"]])) {
        l[["gffVersion"]] <- numeric_version(l[["gffVersion"]])
    }
    ## Parse the head of the GFF file. We can use this to identify ambiguous
    ## sources or genomes with non-standard identifiers.
    lines <- import(
        con = .cacheIt(file),
        format = "lines",
        nMax = 1000L,
        quiet = TRUE
    )
    lines <- lines[!grepl(pattern = "^#", x = lines)]
    df2 <- import(
        con = textConnection(lines),
        format = "tsv",
        colnames = FALSE,
        quiet = TRUE
    )
    if (!isString(l[["provider"]])) {
        if (isSubset(
            x = c(
                "genebuild-last-updated", "genome-build",
                "genome-build-accession", "genome-version"
            ),
            y = df[["key"]]
        )) {
            l[["provider"]] <- "Ensembl"
        } else if (grepl(
            pattern = paste0(
                "^(",
                "ensGene", "|",
                "knownGene", "|",
                "ncbiRefSeq", "|",
                ## Now seeing this in files as of 2023.
                "ncbiRefSeq\\.[0-9]{4}-[0-9]{2}-[0-9]{2}", "|",
                "refGene",
                ")$"
            ),
            x = df2[[2L]][[1L]]
        )) {
            l[["provider"]] <- "UCSC"
        } else if (identical("FlyBase", df2[[2L]][[1L]])) {
            ## FIXME This returns incorrectly for:
            ## https://ftp.ensembl.org/pub/release-110/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.46.110.gtf.gz
            l[["provider"]] <- "FlyBase"
        } else if (identical("WormBase", df2[[2L]][[1L]])) {
            ## FIXME Need to check C. elegans Ensembl file.
            l[["provider"]] <- "WormBase"
        } else if (
            identical(basename(file), "ref-transcripts.gtf") &&
                identical("chrI", df2[[1L]][[1L]])
        ) {
            ## bcbio-nextgen genome, processed with gffutils.
            ## https://github.com/daler/gffutils
            l[["provider"]] <- "UCSC"
        } else {
            abort(sprintf(
                "Failed to detect provider (e.g. {.val %s}) from {.file %s}.",
                "Ensembl", basename(file)
            ))
        }
    }
    assert(isString(l[["provider"]]))
    pattern <- .gffPatterns[[tolower(l[["provider"]])]]
    assert(isString(pattern))
    switch(
        EXPR = l[["provider"]],
        "Ensembl" = {
            if (isTRUE(grepl(pattern = pattern, x = basename(file)))) {
                x <- strMatch(
                    x = basename(file),
                    pattern = pattern
                )[1L, , drop = TRUE]
                if (!isOrganism(l[["organism"]])) {
                    l[["organism"]] <- gsub("_", " ", x[[3L]])
                }
                if (!isString(l[["genomeBuild"]])) {
                    l[["genomeBuild"]] <- x[[4L]]
                }
                if (!isInt(l[["release"]])) {
                    l[["release"]] <- as.integer(x[[5L]])
                }
            }
        },
        "FlyBase" = {
            if (isTRUE(grepl(pattern = pattern, x = basename(file)))) {
                x <- strMatch(
                    x = basename(file),
                    pattern = pattern
                )[1L, , drop = TRUE]
                if (
                    !isString(l[["organism"]]) &&
                    identical(x[[3L]], "dmel")
                ) {
                    l[["organism"]] <- "Drosophila melanogaster"
                }
                if (!isString(l[["release"]])) {
                    l[["release"]] <- x[[5L]]
                }
                if (!isString(l[["genomeBuild"]])) {
                    l[["genomeBuild"]] <- l[["release"]]
                }
            }
        },
        "GENCODE" = {
            if (isTRUE(grepl(pattern = pattern, x = basename(file)))) {
                x <- strMatch(
                    x = basename(file),
                    pattern = pattern
                )[1L, , drop = TRUE]
                if (!isScalar(l[["release"]])) {
                    l[["release"]] <- x[[3L]]
                    if (grepl(
                        pattern = "^[0-9]+",
                        x = l[["release"]]
                    )) {
                        l[["release"]] <-
                            as.integer(l[["release"]])
                    }
                }
            }
        },
        "RefSeq" = {
            if (!isScalar(l[["release"]])) {
                ## e.g. "109.20190125".
                l[["release"]] <- strMatch(
                    x = df[
                        df[["key"]] == "annotation-source",
                        "value",
                        drop = TRUE
                    ],
                    pattern = "^NCBI.+Annotation\\sRelease\\s([.0-9]+)$"
                )[1L, 2L]
            }
        },
        "UCSC" = {
            if (!isString(l[["genomeBuild"]])) {
                l[["genomeBuild"]] <- strMatch(
                    x = basename(file),
                    pattern = pattern
                )[1L, 3L]
            }
        },
        "WormBase" = {
            if (!isString(l[["release"]])) {
                l[["release"]] <- df[
                    df[["key"]] == "genebuild-version",
                    "value",
                    drop = TRUE
                ]
            }
            if (isTRUE(grepl(pattern = pattern, x = basename(file)))) {
                x <- strMatch(
                    x = basename(file),
                    pattern = pattern
                )[1L, , drop = TRUE]
                if (
                    !isString(l[["organism"]]) &&
                    identical(x[[3L]], "c_elegans")
                ) {
                    l[["organism"]] <- "Caenorhabditis elegans"
                }
                if (!isString(l[["release"]])) {
                    l[["release"]] <- x[[5L]]
                }
            }
            if (!isString(l[["genomeBuild"]])) {
                l[["genomeBuild"]] <- l[["release"]]
            }
        }
    )
    ## Attempt to detect the organism from the genome build, if necessary.
    if (!isOrganism(l[["organism"]]) && isString(l[["genomeBuild"]])) {
        l[["organism"]] <- tryCatch(
            expr = detectOrganism(l[["genomeBuild"]]),
            error = function(e) {
                NULL
            }
        )
    }
    ## Attempt to detect the organism from gene identifiers, if necessary.
    if (!isOrganism(l[["organism"]])) {
        match <- strMatch(
            x = lines,
            pattern = "(\t|\\s)gene_id\\s\"([^\"]+)\""
        )
        genes <- unique(na.omit(match[, 3L, drop = TRUE]))
        l[["organism"]] <- tryCatch(
            expr = detectOrganism(genes),
            error = function(e) {
                NULL
            }
        )
    }
    ## Genome-specific organism workarounds.
    if (
        !isOrganism(l[["organism"]]) &&
            any(grepl(pattern = "gene_source \"sgd\"", x = lines))
    ) {
        l[["organism"]] <- "Saccharomyces cerevisiae"
    }
    l <- Filter(f = hasLength, x = l)
    l <- l[sort(names(l))]
    assert(
        isString(l[["file"]]),
        isString(l[["format"]]),
        isString(l[["md5"]]),
        isOrganism(l[["organism"]]),
        isString(l[["provider"]]),
        isString(l[["sha256"]]),
        msg = sprintf("Unsupported GFF: {.file %s}.", basename(file))
    )
    l
}



## Updated 2023-12-04.
.gffGenomeBuild <- function(df) {
    assert(is(df, "DFrame"))
    ## GENCODE files have a description key that contains the genome build.
    if (
        isSubset(c("description", "provider"), df[["key"]]) &&
            identical(
                x = df[which(df[["key"]] == "provider"), "value"],
                y = "GENCODE"
            )
    ) {
        string <- df[df[["key"]] == "description", "value", drop = TRUE]
        ## This edge case applies to GRCh37, lifted over from GRCh38.
        match <- strMatch(x = string, pattern = "mapped to ([^ ]+)")
        if (identical(dim(match), c(1L, 2L)) && isString(match[1L, 2L])) {
            x <- match[1L, 2L]
            return(x)
        }
        match <- strMatch(x = string, pattern = "genome \\(([^\\)]+)\\)")
        assert(
            identical(dim(match), c(1L, 2L)),
            isString(match[1L, 2L])
        )
        x <- match[1L, 2L]
        return(x)
    }
    ## Otherwise we can parse for standard "genome-build" key, which is
    ## supported by Ensembl and RefSeq.
    .getValue <- function(key, df) {
        x <- df[match(x = key, table = df[["key"]]), "value", drop = TRUE]
        if (is.na(x)) {
            return(NULL)
        }
        x
    }
    x <- .getValue(key = "genome-build", df = df)
    if (isString(x)) {
        if (isTRUE(grepl(pattern = " ", x = x))) {
            x <- strsplit(x = x, split = " ", fixed = TRUE)[[1L]]
            x <- x[length(x)]
        }
        return(x)
    }
    NULL
}



## Updated 2021-01-21.
.gffFormat <- function(file) {
    ifelse(
        test = grepl(
            pattern = "^gtf",
            x = fileExt(file),
            ignore.case = TRUE
        ),
        yes = "GTF",
        no = "GFF3"
    )
}



## Updated 2021-01-22.
.gffProvider <- function(df) {
    assert(is(df, "DFrame"))
    annoSource <-
        df[df[["key"]] == "annotation-source", "value", drop = TRUE]
    geneBuildVersion <-
        df[df[["key"]] == "genebuild-version", "value", drop = TRUE]
    provider <-
        df[df[["key"]] == "provider", "value", drop = TRUE]
    if (identical(provider, "GENCODE")) {
        return("GENCODE")
    } else if (isTRUE(grepl(pattern = "^NCBI", x = annoSource))) {
        return("RefSeq")
    } else if (isSubset(
        x = c(
            "genebuild-last-updated",
            "genome-build",
            "genome-build-accession",
            "genome-date",
            "genome-version"
        ),
        y = df[["key"]]
    )) {
        return("Ensembl")
    } else if (isTRUE(grepl(pattern = "^WS[0-9]+$", x = geneBuildVersion))) {
        return("WormBase")
    }
    NULL
}



#' Is the input GFF file supported in the package?
#'
## See `.gffPatterns` for pattern matching details.
#'
#' @note Updated 2023-11-29.
#' @noRd
.isSupportedGff <- function(file) {
    ok <- isString(file)
    if (!ok) {
        return(FALSE)
    }
    denylist <- c(
        "flybase_gff" = paste0(
            "^([a-z0-9]+_)?",
            "^([^-]+)",
            "-([^-]+)",
            "-(r[0-9]+\\.[0-9]+)",
            "\\.gff",
            "(\\.gz)?$"
        ),
        "wormbase_gff" = paste0(
            "^([a-z0-9]+_)?",
            "^([a-z]_[a-z]+)",
            "\\.([A-Z0-9]+)",
            "\\.(WS[0-9]+)",
            "\\.([a-z_]+)",
            "\\.gff3",
            "(\\.gz)?$"
        )
    )
    if (grepl(
        pattern = denylist[["flybase_gff"]],
        x = basename(file)
    )) {
        alertWarning("Use FlyBase GTF instead of GFF.")
        return(FALSE)
    } else if (grepl(
        pattern = denylist[["wormbase_gff"]],
        x = basename(file)
    )) {
        alertWarning("Use WormBase GTF instead of GFF.")
        return(FALSE)
    }
    TRUE
}
