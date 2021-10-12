#' Get metadata about a GFF file
#'
#' @note Updated 2021-08-06.
#' @export
#'
#' @inheritParams AcidRoxygen::params
#'
#' @return `list`.
#'   Containing values, when possible:
#'   - `directives`
#'   - `format`
#'   - `genomeBuild`
#'   - `gffVersion`
#'   - `md5`
#'   - `organism`
#'   - `provider`
#'   - `release`
#'   - `sha256`
#'
#' @seealso
#' - [getGFFDirectives()].
#'
#' @examples
#' url <- pasteURL(
#'     "ftp.ensembl.org",
#'     "pub",
#'     "release-102",
#'     "gtf",
#'     "homo_sapiens",
#'     "Homo_sapiens.GRCh38.102.gtf.gz",
#'     protocol = "ftp"
#' )
#' x <- getGFFMetadata(url)
#' print(x)
getGFFMetadata <- function(file) {
    assert(.isSupportedGFF(file))
    l <- list()
    if (isAFile(file)) {
        file <- realpath(file)
    }
    alert(sprintf("Getting GFF metadata for {.file %s}.", basename(file)))
    l[["file"]] <- file
    l[["md5"]] <- .md5(file)
    l[["sha256"]] <- .sha256(file)
    ## Attempt to get genome build and provider from GFF directives.
    df <- getGFFDirectives(file)
    if (is(df, "DataFrame")) {
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
    if (!isString(l[["provider"]])) {
        lines <- import(
            file = .cacheIt(file),
            format = "lines",
            nMax = 1000L,
            quiet = TRUE
        )
        lines <- lines[!grepl(pattern = "^#", x = lines)]
        lines <- head(lines, n = 10L)
        if (any(grepl(
            pattern = "\t(ensGene|knownGene|ncbiRefSeq|refGene)\t",
            x = lines
        ))) {
            l[["provider"]] <- "UCSC"
        } else if (any(grepl(
            pattern = "\tFlyBase\t", x = lines
        ))) {
            l[["provider"]] <- "FlyBase"
        } else if (any(grepl(
            pattern = "\tWormBase\t", x = lines
        ))) {
            l[["provider"]] <- "WormBase"
        } else if (any(grepl(
            ## NOTE This will currently miss genomes with gene identifiers that
            ## aren't prefixed with "ENS".
            pattern = "\tgene_id \"ENS.*G[0-9]{11}", x = lines
        ))) {
            l[["provider"]] <- "Ensembl"
        } else {
            abort(sprintf(
                "Failed to detect provider (e.g. {.val %s}) from {.file %s}.",
                "Ensembl", basename(file)
            ))
        }
    }
    ## Attempt to parse file names for useful values.
    if (isString(l[["provider"]])) {
        pattern <- .gffPatterns[[tolower(l[["provider"]])]]
        switch(
            EXPR = l[["provider"]],
            "Ensembl" = {
                if (isTRUE(grepl(pattern = pattern, x = basename(file)))) {
                    x <- str_match(
                        string = basename(file),
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
                    x <- str_match(
                        string = basename(file),
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
                    x <- str_match(
                        string = basename(file),
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
                    l[["release"]] <- str_match(
                        string = df[
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
                    l[["genomeBuild"]] <- str_match(
                        string = basename(file),
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
                    x <- str_match(
                        string = basename(file),
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
    }
    ## Attempt to detect the organism from the genome build, if necessary.
    if (
        !isOrganism(l[["organism"]]) &&
        isString(l[["genomeBuild"]])
    ) {
        l[["organism"]] <- tryCatch(
            expr = detectOrganism(l[["genomeBuild"]]),
            error = function(e) NULL
        )
    }
    ## Attempt to detect the organism from gene identifiers, if necessary.
    if (!isOrganism(l[["organism"]])) {
        match <- str_match(
            string = lines,
            pattern = "\tgene_id \"([^\"]+)\""
        )
        genes <- na.omit(match[, 2L, drop = TRUE])
        l[["organism"]] <- tryCatch(
            expr = detectOrganism(genes),
            error = function(e) NULL
        )
    }
    l <- Filter(f = hasLength, x = l)
    l <- l[sort(names(l))]
    assert(
        isString(l[["file"]]),
        isString(l[["format"]]),
        isString(l[["md5"]]),
        isOrganism(l[["organism"]]),
        isString(l[["provider"]]),
        isString(l[["sha256"]])
    )
    l
}



## Updated 2021-08-05.
.gffGenomeBuild <- function(df) {
    assert(is(df, "DataFrame"))
    ## GENCODE files have a description key that contains the genome build.
    if (
        identical(
            x = df[which(df[["key"]] == "provider"), "value"],
            y = "GENCODE"
        ) &&
        isTRUE("description" %in% df[["key"]])
    ) {
        string <- df[df[["key"]] == "description", "value", drop = TRUE]
        x <- str_match(
            string = string,
            pattern = "genome \\(([^\\)]+)\\)"
        )[1L, 2L]
        if (isString(x)) {
            return(x)
        }
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
    assert(is(df, "DataFrame"))
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
