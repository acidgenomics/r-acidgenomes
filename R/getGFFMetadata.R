#' Get metadata about a GFF file
#'
#' @note Updated 2021-01-22.
#' @export
#'
#' @inheritParams AcidRoxygen::params
#'
#' @return `list`.
#'   Containing values, if defined:
#'   - `directives`
#'   - `format`
#'   - `genomeBuild`
#'   - `gffVersion`
#'   - `organism`
#'   - `provider`
#'   - `providerVersion`
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
    l <- list()
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
    if (!isString(l[["provider"]])) {
        lines <- import(
            file = .cacheIt(file),
            format = "lines",
            nMax = 5L,
            quiet = TRUE
        )
        if (any(grepl(
            pattern = "\t(ensGene|knownGene|ncbiRefSeq|refGene)\t",
            x = lines
        ))) {
            l[["provider"]] <- "UCSC"
        } else if (any(grepl(pattern = "\tFlyBase\t", x = lines))) {
            l[["provider"]] <- "FlyBase"
        }
    }
    if (!isString(l[["format"]])) {
        l[["format"]] <- .gffFormat(file)
    }
    assert(isSubset(l[["format"]], .gffFormats))
    if (isString(l[["gffVersion"]])) {
        l[["gffVersion"]] <- numeric_version(l[["gffVersion"]])
    }
    ## Attempt to parse file names for useful values.
    if (isString(l[["provider"]])) {
        switch(
            EXPR = l[["provider"]],
            "Ensembl" = {
                pattern <- paste0(
                    "^([a-z0-9]+_)?",                # BiocFileCache
                    "^([A-Z][a-z]+_[a-z]+)",         # "Homo_sapiens"
                    "\\.([A-Za-z0-9]+)",             # "GRCh38"
                    "\\.([0-9]+)",                   # "102"
                    "(\\.chr_patch_hapl_scaff)?",
                    "\\.(gff3|gtf)",
                    "(\\.gz)?$"
                )
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
                    if (!isInt(l[["providerVersion"]])) {
                        l[["providerVersion"]] <- as.integer(x[[5L]])
                    }
                }
            },
            "FlyBase" = {
                pattern <- paste0(
                    "^([a-z0-9]+_)?",       # BiocFileCache
                    "^([^-]+)",             # "dmel"
                    "-([^-]+)",             # "all"
                    "-(r[0-9]+\\.[0-9]+)",  # "r6.37"
                    "\\.(gff|gtf)",
                    "(\\.gz)?$"
                )
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
                }
                if (!isString(l[["providerVersion"]])) {
                    l[["providerVersion"]] <- x[[5L]]
                }
                if (!isString(l[["genomeBuild"]])) {
                    ## This matches the convention defined in GFF.
                    l[["genomeBuild"]] <- paste(
                        provider, providerVersion
                    )
                }
            },
            "GENCODE" = {
                pattern <- paste0(
                    "^([a-z0-9]+_)?",                # BiocFileCache
                    "gencode",
                    "\\.v([M0-9]+)",                 # 36 (human) / M25 (mouse)
                    "(lift37)?",                     # GRCh37-specific
                    "\\.annotation",
                    "\\.(gff3|gtf)",
                    "(\\.gz)?$"
                )
                if (isTRUE(grepl(pattern = pattern, x = basename(file)))) {
                    x <- str_match(
                        string = basename(file),
                        pattern = pattern
                    )[1L, , drop = TRUE]
                    if (!isScalar(l[["providerVersion"]])) {
                        l[["providerVersion"]] <- x[[3L]]
                        if (grepl(
                            pattern = "^[0-9]+",
                            x = l[["providerVersion"]]
                        )) {
                            l[["providerVersion"]] <-
                                as.integer(l[["providerVersion"]])
                        }
                    }
                }
            },
            "RefSeq" = {
                if (!isInt(l[["providerVersion"]])) {
                    l[["providerVersion"]] <- as.integer(str_match(
                        string = df[
                            df[["key"]] == "annotation-source",
                            "value",
                            drop = TRUE
                        ],
                        pattern = "^NCBI.+Annotation\\sRelease\\s([0-9]+)$"
                    )[1L, 2L])
                }
            },
            "UCSC" = {
                if (!isString(l[["genomeBuild"]])) {
                    l[["genomeBuild"]] <- str_match(
                        string = basename(file),
                        pattern = paste0(
                            "^([0-9a-z]_)?",            # BiocFileCache.
                            "([a-z]+[A-Za-z]+[0-9]+)",  # "hg38"
                            "\\."
                        )
                    )[1L, 3L]
                }
            }
        )
    }
    ## Attempt to get the organism from the genome build, if necessary.
    if (!isOrganism(l[["organism"]])) {
        l[["organism"]] <- tryCatch(
            expr = detectOrganism(l[["genomeBuild"]]),
            error = function(e) NULL
        )
    }
    l <- Filter(f = hasLength, x = l)
    l <- l[sort(names(l))]
    l
}



## Updated 2021-01-20.
.gffGenomeBuild <- function(df) {
    assert(is(df, "DataFrame"))
    ## GENCODE files have a description key that contains the genome build.
    if (isTRUE("description" %in% df[["key"]])) {
        string <- df[df[["key"]] == "description", "value", drop = TRUE]
        x <- str_match(
            string = string,
            pattern = "genome \\(([^\\)]+)\\)"
        )[1L, 2L]
        if (isString(x)) return(x)
    }
    ## Otherwise we can parse for standard "genome-build" key, which is
    ## supported by Ensembl and RefSeq.
    .getValue <- function(key) {
        x <- df[match(x = key, table = df[["key"]]), "value", drop = TRUE]
        if (is.na(x)) return(NULL)
        x
    }
    x <- .getValue("genome-build")
    if (isString(x)) return(x)
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
    }
    NULL
}
