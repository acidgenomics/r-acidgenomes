## Updated 2021-01-22.
.gffMetadata <- function(file) {
    l <- list()
    l[["type"]] <- .gffType(file)
    ## Attempt to get genome build and source from GFF directives.
    df <- getGFFDirectives(file, nMax = 2000L)
    if (is(df, "DataFrame")) {
        l[["genomeBuild"]] <- .gffGenomeBuild(df)
        l[["source"]] <- .gffSource(df)
    } else {
        ## Otherwise, fall back to attempting to get the source from the
        ## first lines of the file.
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
            l[["source"]] <- "UCSC"
        }
    }
    ## Attempt to parse file names for useful values.
    if (!is.null(l[["source"]])) {
        switch(
            EXPR = l[["source"]],
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
                    if (is.null(l[["organism"]])) {
                        l[["organism"]] <- gsub("_", " ", x[[3L]])
                    }
                    if (is.null(l[["genomeBuild"]])) {
                        l[["genomeBuild"]] <- x[[4L]]
                    }
                    if (is.null(l[["release"]])) {
                        l[["release"]] <- as.integer(x[[5L]])
                    }
                }
            },
            "UCSC" = {
                if (is.null(l[["genomeBuild"]])) {
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
    if (is.null(l[["organism"]])) {
        l[["organism"]] <- tryCatch(
            expr = detectOrganism(l[["genomeBuild"]]),
            error = function(e) NULL
        )
    }
    l <- Filter(f = Negate(is.null), x = l)
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



## Updated 2021-01-20.
.gffSource <- function(df) {
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



## Updated 2021-01-21.
.gffType <- function(file) {
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
