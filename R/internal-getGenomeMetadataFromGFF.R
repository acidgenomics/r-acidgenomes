#' Get genome metadata from a GFF file
#'
#' @note Updated 2021-01-18.
#' @noRd
#'
#' @examples
#' files <- c(
#'     pasteURL(
#'         "ftp.ensembl.org",
#'         "pub",
#'         "release-102",
#'         "gtf",
#'         "homo_sapiens",
#'         "Homo_sapiens.GRCh38.102.gtf.gz",
#'         protocol = "ftp"
#'     ),
#'     pasteURL(
#'         "ftp.ensembl.org",
#'         "pub",
#'         "release-102",
#'         "gff3",
#'         "homo_sapiens",
#'         "Homo_sapiens.GRCh38.102.gff3.gz",
#'         protocol = "ftp"
#'     ),
#'     pasteURL(
#'         "ftp.ensembl.org",
#'         "pub",
#'         "grch37",
#'         "release-102",
#'         "gtf",
#'         "homo_sapiens",
#'         "Homo_sapiens.GRCh37.87.gtf.gz",
#'         protocol = "ftp"
#'     ),
#'     pasteURL(
#'         "ftp.ensembl.org",
#'         "pub",
#'         "grch37",
#'         "release-102",
#'         "gff3",
#'         "homo_sapiens",
#'         "Homo_sapiens.GRCh37.87.gff3.gz",
#'         protocol = "ftp"
#'     ),
#'     pasteURL(
#'         "ftp.ebi.ac.uk",
#'         "pub",
#'         "databases",
#'         "gencode",
#'         "Gencode_human",
#'         "release_36",
#'         "gencode.v36.annotation.gtf.gz",
#'         protocol = "ftp"
#'     ),
#'     pasteURL(
#'         "ftp.ebi.ac.uk",
#'         "pub",
#'         "databases",
#'         "gencode",
#'         "Gencode_human",
#'         "release_36",
#'         "gencode.v36.annotation.gff3.gz",
#'         protocol = "ftp"
#'     ),
#'     pasteURL(
#'         "ftp.ebi.ac.uk",
#'         "pub",
#'         "databases",
#'         "gencode",
#'         "Gencode_human",
#'         "release_36",
#'         "GRCh37_mapping",
#'         "gencode.v36lift37.annotation.gtf.gz",
#'         protocol = "ftp"
#'     ),
#'     pasteURL(
#'         "ftp.ebi.ac.uk",
#'         "pub",
#'         "databases",
#'         "gencode",
#'         "Gencode_human",
#'         "release_36",
#'         "GRCh37_mapping",
#'         "gencode.v36lift37.annotation.gff3.gz",
#'         protocol = "ftp"
#'     ),
#'     pasteURL(
#'         "ftp.ebi.ac.uk",
#'         "pub",
#'         "databases",
#'         "gencode",
#'         "Gencode_mouse",
#'         "release_M25",
#'         "gencode.vM25.annotation.gtf.gz",
#'         protocol = "ftp"
#'     ),
#'     pasteURL(
#'         "ftp.ebi.ac.uk",
#'         "pub",
#'         "databases",
#'         "gencode",
#'         "Gencode_mouse",
#'         "release_M25",
#'         "gencode.vM25.annotation.gff3.gz",
#'         protocol = "ftp"
#'     ),
#'     pasteURL(
#'         "ftp.ncbi.nlm.nih.gov",
#'         "genomes",
#'         "refseq",
#'         "vertebrate_mammalian",
#'         "Homo_sapiens",
#'         "all_assembly_versions",
#'         "GCF_000001405.38_GRCh38.p12",
#'         "GCF_000001405.38_GRCh38.p12_genomic.gtf.gz",
#'         protocol = "ftp"
#'     ),
#'     pasteURL(
#'         "ftp.ncbi.nlm.nih.gov",
#'         "genomes",
#'         "refseq",
#'         "vertebrate_mammalian",
#'         "Homo_sapiens",
#'         "all_assembly_versions",
#'         "GCF_000001405.38_GRCh38.p12",
#'         "GCF_000001405.38_GRCh38.p12_genomic.gff.gz",
#'         protocol = "ftp"
#'     ),
#'     ## Note that this file doesn't contain any metadata comments.
#'     pasteURL(
#'         "ftp.flybase.net",
#'         "releases",
#'         "FB2020_06",
#'         "dmel_r6.37",
#'         "gtf",
#'         "dmel-all-r6.37.gtf.gz",
#'         protocol = "ftp"
#'     ),
#'     ## This file is very large and slow to parse.
#'     pasteURL(
#'         "ftp.flybase.net",
#'         "releases",
#'         "FB2020_06",
#'         "dmel_r6.37",
#'         "gff",
#'         "dmel-all-r6.37.gff.gz",
#'         protocol = "ftp"
#'     ),
#'     pasteURL(
#'         "ftp.wormbase.org",
#'         "pub",
#'         "wormbase",
#'         "releases",
#'         "WS279",
#'         "species",
#'         "c_elegans",
#'         "PRJNA13758",
#'         "c_elegans.PRJNA13758.WS279.canonical_geneset.gtf.gz",
#'         protocol = "ftp"
#'     ),
#'     pasteURL(
#'         "ftp.wormbase.org",
#'         "pub",
#'         "wormbase",
#'         "releases",
#'         "WS279",
#'         "species",
#'         "c_elegans",
#'         "PRJNA13758",
#'         "c_elegans.PRJNA13758.WS279.annotations.gff3.gz",
#'         protocol = "ftp"
#'     )
#' )
#' lapply(X = files, FUN = .getGenomeMetadataFromGFF)

.getGenomeMetadataFromGFF <- function(file) {
    meta <- .getGFFMetadata(file, nMax = 100L)
    assert(is(meta, "DataFrame"))

    ## FIXME DROP SEQUENCE REGIONS
    ## FIXME SORT ALPHABETICALLY AND MAKE UNIQUE?
    ## FIXME ENSURE KEYS ARE UNIQUE...

    ## FIXME How to handle sequence region comments here?
    patterns <- c(
        "ensembl" = paste0(
            "^([a-z0-9]+_)?",          # temp prefix from BiocFileCache.
            "([A-Z][a-z]+_[a-z]+)",    # organism (e.g. "Homo_sapiens").
            "\\.([A-Za-z0-9]+)",       # genomeVersion (e.g. "GRCh38").
            "\\.([0-9]+)",             # (Ensembl release) version (e.g. "102").
            "\\.(gff3|gtf)",
            "(\\.gz)?$"
        ),
        "gencode" = paste0(
            "^([a-z0-9]+_)?",          # temp prefix from BiocFileCache.
            "gencode",
            "\\.v([M0-9]+)",           # v32 (human) or vM25 (mouse).
            "(lift37)?",               # GRCh37 liftover.
            "\\.annotation",
            "\\.(gff3|gtf)",
            "(\\.gz)?$"
        ),
        "refseq" = paste0(
            "^([a-z0-9]+_)?",          # temp prefix from BiocFileCache.
            "(GCF_[0-9]+\\.[0-9]+)",   # Accession (e.g. GCF_000001405.38).
            "_(.+)",                   # Genome build (e.g. GRCh38.p12).
            "_genomic",
            "\\.(gff|gtf)",
            "(\\.gz)?$"
        )
    )
    ## Loop across the genome source patterns, and see if we get a hit.
    hits <- bapply(
        X = patterns,
        string = basename(file),
        FUN = function(pattern, string) {
            isTRUE(grepl(pattern = pattern, x = string))
        }
    )
    if (!any(hits)) {
        ## nocov start
        stop(sprintf(
            paste(
                "Failed to match genome metadata from file name ('%s').",
                "Define values manually: %s."
            ),
            basename(file),
            toString(c("organism", "genomeBuild", "release"))
        ))
        ## nocov end
    }
    source <- names(patterns)[hits][[1L]]
    match <- str_match(
        string = basename(file),
        pattern = patterns[[source]]
    )
    match <- match[1L, , drop = TRUE]
    switch(
        EXPR = source,
        "ensembl" = {
            alert("Detecting Ensembl genome metadata from file name.")
            organism <- gsub(
                pattern = "_",
                replacement = " ",
                x = match[[3L]]
            )
            genomeBuild <- match[[4L]]
            release <- as.integer(match[[5L]])
        },
        "gencode" = {
            alert("Detecting GENCODE genome metadata from file name.")
            release <- match[[3L]]
            if (grepl("^M", release)) {
                organism <- "Mus musculus"
            } else {
                organism <- "Homo sapiens"
                release <- as.integer(release)
                ## GRCh38
                if (identical(match[[4L]], "lift37")) {
                    genomeBuild <- "GRCh37"
                }
            }
            ## Assembly (genome build) is documented in the first commented
            ## lines of the file (line 1 for GTF; line 2 for GFF3).
            if (is.null(genomeBuild)) {
                lines <- import(
                    file = file,
                    format = "lines",
                    nMax = 2L,
                    quiet = TRUE
                )
                genomeBuild <- str_match(
                    string = grep(
                        ## FIXME TIGHTEN UP THIS MATCH.
                        pattern = "description:",
                        x = lines,
                        value = TRUE
                    ),
                    pattern = " genome \\(([^\\)]+)\\),"
                )[1L, 2L]
            }
        },
        "refseq" = {
            lines <- import(
                file = file,
                format = "lines",
                nMax = 4L,
                quiet = TRUE
            )
            release <- as.integer(str_match(
                string = grep(
                    pattern = "^#!annotation-source",
                    x = lines,
                    value = TRUE
                ),
                pattern = "Annotation Release ([0-9]+)"
            )[1L, 2L])
            genomeBuild <- match[[4L]]
            release <- NA_integer_
        }
    )
    assert(
        isOrganism(organism),
        isString(genomeBuild),
        isInt(release) || isString(release)
    )
    list(
        "comments" = lines,
        "organism" = organism,
        "genomeBuild" = genomeBuild,
        "accession" = "FIXME",  # GCA_000001405.28
        "release" = release,
        "source" = "FIXME"
    )
}
