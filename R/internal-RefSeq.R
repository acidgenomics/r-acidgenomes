## nolint start

#' Get the RefSeq assembly metadata
#'
#' @section Stable reference assembly:
#' The latest reference assembly, linked under the `reference/`
#' subdirectory changes over time and is not considered "stable". For improved
#' reproduciblity, track a reference one version behind
#' (e.g. use "GCF_000001405.39_GRCh38.p12" instead of current
#' "GCF_000001405.39_GRCh38.p13" build).
#'
#' This approach to maintaining a current reference build differs on NCBI RefSeq
#' compared to other sources such as Ensembl and GENCODE, which track a stable
#' release as the latest "reference" assembly.
#'
#' @note Updated 2021-01-25.
#' @noRd
#'
#' @param file `character(1)`.
#'   RefSeq assembly summary file or URL.
#'
#' @seealso
#' - [File format details](ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt).
#'
#' @return Named `character`.
#'
#' @examples
#' file <- pasteURL(
#'     "ftp.ncbi.nlm.nih.gov",
#'     "genomes",
#'     "refseq",
#'     "vertebrate_mammalian",
#'     "Homo_sapiens",
#'     "assembly_summary.txt",
#'     protocol = "ftp"
#' )
#' x <- .getRefSeqAssemblySummary(file)
#' names(x)

## nolint end

.getRefSeqAssemblySummary <-
    function(file) {
        pattern <- "assembly_summary.txt"
        assert(
            isString(file),
            isMatchingFixed(pattern = pattern, x = basename(file))
        )
        file <- .cacheIt(file)
        lines <- import(
            file = file,
            format = "lines",
            skip = 1L,
            quiet = TRUE
        )
        names <- strsplit(
            x = sub(pattern = "^#\\s", replacement = "", x = lines[[1L]]),
            split = "\\t"
        )[[1L]]
        values <- strsplit(x = lines[[2L]], split = "\\t")[[1L]]
        x <- as.character(values[seq_len(20L)])
        names(x) <- names[seq_len(20L)]
        x <- x[nzchar(x)]
        x
    }



#' Get the RefSeq base genome URL for an organism
#'
#' @note Updated 2021-01-08.
#' @noRd
#'
#' @examples
#' .getRefSeqGenomeURL(
#'     organism = "Homo sapiens",
#'     taxonomicGroup = "vertebrate_mammalian"
#' )
.getRefSeqGenomeURL <- function(
    organism,
    taxonomicGroup = NULL,
    quiet = FALSE
) {
    assert(
        isOrganism(organism),
        isString(taxonomicGroup, nullOK = TRUE),
        isFlag(quiet)
    )
    if (isFALSE(quiet)) {
        alert(sprintf(
            "Locating {.emph %s} genome on RefSeq FTP server.",
            organism
        ))
    }
    baseURL <- "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq"
    if (is.null(taxonomicGroup)) {
        alertWarning(sprintf(
            "Set {.var %s} manually to speed up this step.",
            "taxonimicGroup"
        ))
        taxonomicGroups <- getURLDirList(url = baseURL)
        keep <- grepl(pattern = "^[a-z_]+$", x = taxonomicGroups)
        taxonomicGroups <- sort(taxonomicGroups[keep])
        list <- bplapply(
            X = taxonomicGroups,
            baseURL = baseURL,
            FUN = function(taxonomicGroup, baseURL) {
                url <- pasteURL(baseURL, taxonomicGroup)
                x <- getURLDirList(url = url)
                keep <- grepl(pattern = "^[A-Z][a-z]+_[a-z]+$", x = x)
                x <- sort(x[keep])
                x
            }
        )
        names(list) <- taxonomicGroups
        match <- vapply(
            X = list,
            organism = gsub(pattern = " ", replacement = "_", x = organism),
            FUN = function(strings, organism) {
                isSubset(x = organism, y = strings)
            },
            FUN.VALUE = logical(1L),
            USE.NAMES = TRUE
        )
        taxonomicGroup <- names(match)[match]
        assert(isString(taxonomicGroup))
    }
    url <- pasteURL(
        baseURL,
        taxonomicGroup,
        gsub(pattern = " ", replacement = "_", x = organism)
    )
    if (isFALSE(quiet)) {
        dl(c("URL" = url))
    }
    url
}



#' Get RefSeq genome assembly seqinfo
#'
#' Parse the assembly report file to get `seqlengths` per chromosome.
#'
#' @note Updated 2021-01-25.
#' @noRd
#'
#' @param file `character(1)`.
#'   RefSeq GFF file.
#'
#' @return `Seqinfo`.
#'
#' @seealso
#' - `tximeta:::gtf2RefSeq`.
#'
#' @examples
#' file <- pasteURL(
#'     "ftp.ncbi.nlm.nih.gov",
#'     "genomes",
#'     "refseq",
#'     "vertebrate_mammalian",
#'     "Homo_sapiens",
#'     "all_assembly_versions",
#'     "GCF_000001405.38_GRCh38.p12",
#'     "GCF_000001405.38_GRCh38.p12_genomic.gff.gz",
#'     protocol = "ftp"
#' )
#' seqinfo <- .getRefSeqSeqinfo(file)
#' print(seqinfo)
.getRefSeqSeqinfo <- function(file) {
    ## Locate the "*_assembly_report.txt" file from the GFF file path.
    file <- .locateRefSeqAssemblyReport(file)
    file <- .cacheIt(file)
    pattern <- "^([a-z0-9]+_)?GCF_[0-9]+\\.[0-9]+_(.+)_assembly_report\\.txt$"
    assert(isMatchingRegex(pattern = pattern, x = basename(file)))
    ## e.g. GRCh38.p13, which is the format Seqinfo expects.
    ## Refer to GenomeInfoDb documentation for details on NCBI.
    genomeBuild <- sub(
        pattern = pattern,
        replacement = "\\2",
        x = basename(file)
    )
    df <- import(
        file = file,
        format = "tsv",
        colnames = c(
            "sequenceName",
            "sequenceRole",
            "assignedMolecule",
            "assignedMoleculeLocation",
            "genbankAccn",
            "relationship",
            "refseqAccn",
            "assemblyUnit",
            "sequenceLength",
            "ucscStyleName"
        ),
        comment = "#"
    )
    cols <- c("refseqAccn", "sequenceLength")
    df <- df[, cols]
    df <- df[complete.cases(df), ]
    seqnames <- df[["refseqAccn"]]
    seqlengths <- df[["sequenceLength"]]
    assert(
        !any(is.na(seqnames)),
        !any(is.na(seqlengths)),
        hasNoDuplicates(seqnames)
    )
    seq <- Seqinfo(
        seqnames = seqnames,
        seqlengths = seqlengths,
        isCircular = NA,
        genome = genomeBuild
    )
    assert(is(seq, "Seqinfo"))
    validObject(seq)
    seq
}



## FIXME This needs a fix for:
##
## ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/
## GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_assembly_report.txt
##
## from:
##
## ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/
## GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/
## GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gff.gz

#' Locate RefSeq assembly report, from GFF file
#'
#' @note Updated 2021-01-27.
#' @noRd
#'
#' @examples
#' ## RefSeq FTP URL.
#' file <- pasteURL(
#'     "ftp.ncbi.nlm.nih.gov",
#'     "genomes",
#'     "refseq",
#'     "vertebrate_mammalian",
#'     "Homo_sapiens",
#'     "all_assembly_versions",
#'     "GCF_000001405.38_GRCh38.p12",
#'     "GCF_000001405.38_GRCh38.p12_genomic.gff.gz",
#'     protocol = "ftp"
#' )
#' x <- .locateRefSeqAssemblyReport(file)
#' print(x)
#'
#' ## `download-refseq-genome` convention.
#' file <- file.path(
#'     "homo-sapiens-gcf-000001405-39-grch38-p13-refseq-204",
#'     "annotation.gff3.gz"
#' )
#' x <- .locateRefSeqAssemblyReport(file)
#' print(x)
.locateRefSeqAssemblyReport <- function(file) {
    if (isAFile(file)) {
        file <- realpath(file)
    }
    assert(isMatchingRegex(
        pattern = .gffPatterns[["refseq"]],
        x = basename(file)
    ))
    reportBasename <- sub(
        pattern = "_genomic\\.(gff|gtf)(\\.gz)?$",
        replacement = "_assembly_report\\.txt",
        basename(file)
    )
    if (isAURL(file)) {
        x <- pasteURL(dirname(file), reportBasename)
        return(x)
    }
    ## Full local FTP pulldown.
    x <- file.path(dirname(file), reportBasename)
    if (isAFile(x)) return(x)
    ## `download-refseq-genome` download.
    x <- file.path(dirname(dirname(file)), "metadata", reportBasename)
    if (isAFile(x)) return(x)
    stop(sprintf("Failed to locate RefSeq assembly report from '%s'.", file))
}
