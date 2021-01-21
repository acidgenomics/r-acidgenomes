## Updated 2021-01-21.
.gffMetadata <- function(file) {
    genomeBuild <- NULL
    source <- NULL
    type <- .gffType(file)
    df <- getGFFMetadata(file, nMax = 2000L)
    if (is(df, "DataFrame")) {
        genomeBuild <- .gffGenomeBuild(df)
        source <- .gffSource(df)
    } else {
        lines <- import(
            file = .cacheIt(file),
            format = "lines",
            nMax = 5L,
            quiet = TRUE
        )
        ## UCSC GTF file detection.
        if (any(grepl(
            pattern = paste0(
                "\t(",
                "ensGene",
                "|knownGene",
                "|ncbiRefSeq",
                "|refGene",
                ")\t"
            ),
            x = lines
        ))) {
            source <- "UCSC"
            genomeBuild <- str_match(
                string = basename(file),
                pattern = paste0(
                    "^([0-9a-z]_)?",            # BiocFileCache prefix.
                    "([a-z]+[A-Za-z]+[0-9]+)",  # Genome build (e.g. hg38).
                    "\\."
                )
            )[1L, 3L]
        }
    }
    list(
        "genomeBuild" = genomeBuild,
        "source" = source,
        "type" = type
    )
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
