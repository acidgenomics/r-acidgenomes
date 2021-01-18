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
#'     pasteURL(
#'         "ftp.flybase.net",
#'         "releases",
#'         "FB2020_06",
#'         "dmel_r6.37",
#'         "gtf",
#'         "dmel-all-r6.37.gtf.gz",
#'         protocol = "ftp"
#'     ),
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
#'         "c_elegans.PRJNA13758.279.canonical_geneset.gtf.gz",
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
    file <- .cacheIt(file)
    patterns <- c(
        "ensembl" = paste0(
            "^([a-z0-9]+_)?",        # temp prefix from BiocFileCache.
            "([A-Z][a-z]+_[a-z]+)",  # organism (e.g. "Homo_sapiens").
            "\\.([A-Za-z0-9]+)",     # genomeVersion (e.g. "GRCh38").
            "\\.([0-9]+)",           # (Ensembl release) version (e.g. "102").
            "\\.(gff3|gtf)",
            "(\\.gz)?$"
        ),
        "gencode" = paste0(
            "^([a-z0-9]+_)?",        # temp prefix from BiocFileCache.
            "gencode",
            "\\.v([M0-9]+)",         # v32 (human) or vM25 (mouse).
            "(lift37)?",             # GRCh37 liftover.
            "\\.annotation",
            "\\.(gff3|gtf)",
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
    switch(
        EXPR = source,
        "ensembl" = {
            alert("Detecting Ensembl genome metadata from file name.")
            match <- str_match(
                string = basename(file),
                pattern = patterns[["ensembl"]]
            )
            organism <- gsub(
                pattern = "_",
                replacement = " ",
                x = match[1L, 3L]
            )
            genomeBuild <- match[1L, 4L]
            release <- as.integer(match[1L, 5L])
        },
        "gencode" = {
            alert("Detecting GENCODE genome metadata from file name.")
            match <- str_match(
                string = basename(file),
                pattern = patterns[["gencode"]]
            )
            release <- match[1L, 3L]
            if (grepl("^M", release)) {
                organism <- "Mus musculus"
            } else {
                organism <- "Homo sapiens"
                release <- as.integer(release)
                ## GRCh38
                if (identical(match[1L, 4L], "lift37")) {
                    genomeBuild <- "GRCh37"
                }
            }
            ## Assembly (genome build) is documented in the first commented
            ## lines of the file (line 1 for GTF; line 2 for GFF3).
            if (is.null(genomeBuild)) {
                x <- import(
                    file = tmpfile,
                    format = "lines",
                    nMax = 2L,
                    quiet = TRUE
                )
                x <- grep(pattern = "description:", x = x, value = TRUE)
                match <- str_match(
                    string = x,
                    pattern = "\\sgenome\\s\\(([^\\)]+)\\),"
                )
                genomeBuild <- match[1L, 2L]
            }
        }
    )
    assert(
        isOrganism(organism),
        isString(genomeBuild),
        isInt(release) || isString(release)
    )
    list(
        "organism" = organism,
        "genomeBuild" = genomeBuild,
        "release" = release
    )
}
