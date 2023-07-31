#' Match NCBI taxonomic group for gene info or RefSeq.
#'
#' @note Updated 2022-09-22.
#' @noRd
.matchNcbiTaxonomicGroup <-
    function(organism,
             mode = c("geneInfo", "refseq")) {
        assert(isOrganism(organism))
        mode <- match.arg(mode)
        baseUrl <- switch(
            EXPR = mode,
            "geneInfo" = pasteURL(
                "ftp.ncbi.nih.gov",
                "gene", "DATA", "GENE_INFO",
                protocol = "ftp"
            ),
            "refseq" = pasteURL(
                "ftp.ncbi.nlm.nih.gov",
                "genomes", "refseq",
                protocol = "ftp"
            )
        )
        if (isSubset(
            x = organism,
            y = c(
                "Danio rerio", # Zebrafish
                "Homo sapiens", # Human
                "Mus musculus", # Mouse
                "Rattus norvegicus" # Rat
            )
        )) {
            return(switch(
                EXPR = mode,
                "geneInfo" = "Mammalia",
                "refseq" = "vertebrate_mammalian"
            ))
        } else if (isSubset(
            x = organism,
            y = c(
                "Caenorhabditis elegans", # Worm
                "Drosophila melanogaster" # Fruitfly
            )
        )) {
            return(switch(
                EXPR = mode,
                "geneInfo" = "Invertebrates",
                "refseq" = "invertebrate"
            ))
        } else if (isSubset(
            x = organism,
            y = "Saccharomyces cerevisiae" # Yeast
        )) {
            return(switch(
                EXPR = mode,
                "geneInfo" = "Fungi",
                "refseq" = "fungi"
            ))
        }
        alertWarning(sprintf(
            paste(
                "Detecting taxonomic group from {.var %s} at {.url %s}.",
                "Set {.var %s} manually to speed up this step."
            ),
            "organism", baseUrl, "taxonimicGroup"
        ))
        x <- getURLDirList(url = baseUrl)
        pattern <- switch(
            EXPR = mode,
            "geneInfo" = "^[A-Z][A-Za-z_-]+$",
            "refseq" = "^[a-z_]+$"
        )
        keep <- grepl(pattern = pattern, x = x)
        groups <- sort(x[keep])
        list <- lapply(
            X = groups,
            baseUrl = baseUrl,
            mode = mode,
            FUN = function(group, baseUrl, mode) {
                url <- pasteURL(baseUrl, group)
                x <- getURLDirList(url = url)
                switch(
                    EXPR = mode,
                    "geneInfo" = {
                        keep <- grepl(
                            pattern = "^[A-Z][a-z]+_[a-z]+\\.gene_info\\.gz$",
                            x = x
                        )
                        x <- x[keep]
                        x <- gsub(
                            pattern = "\\.gene_info\\.gz$",
                            replacement = "",
                            x = x
                        )
                    },
                    "refseq" = {
                        keep <- grepl(
                            pattern = "^[A-Z][a-z]+_[a-z]+$",
                            x = x
                        )
                        x <- x[keep]
                    }
                )
                x <- sort(x)
                x
            }
        )
        names(list) <- groups
        match <- vapply(
            X = list,
            organism = gsub(pattern = " ", replacement = "_", x = organism),
            FUN = function(strings, organism) {
                isSubset(x = organism, y = strings)
            },
            FUN.VALUE = logical(1L),
            USE.NAMES = TRUE
        )
        out <- names(match)[match]
        assert(isString(out))
        out
    }



## nolint start

#' Get the RefSeq assembly metadata
#'
#' @section Stable reference assembly:
#'
#' The latest reference assembly, linked under the `reference/` subdirectory
#' changes over time and is not considered "stable". For improved
#' reproduciblity, track a reference one version behind
#' (e.g. use "GCF_000001405.39_GRCh38.p12" instead of current
#' "GCF_000001405.39_GRCh38.p13" build).
#'
#' This approach to maintaining a current reference build differs on NCBI RefSeq
#' compared to other sources such as Ensembl and GENCODE, which track a stable
#' release as the latest "reference" assembly.
#'
#' @note Updated 2021-02-12.
#' @noRd
#'
#' @param file `character(1)`.
#' RefSeq assembly summary file or URL.
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
            isMatchingFixed(x = basename(file), pattern = pattern)
        )
        file <- .cacheIt(file)
        lines <- import(
            con = file,
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
#' @note Updated 2023-04-14.
#' @noRd
#'
#' @examples
#' .getRefSeqGenomeUrl(
#'     organism = "Homo sapiens",
#'     taxonomicGroup = "vertebrate_mammalian"
#' )
.getRefSeqGenomeUrl <-
    function(organism,
             taxonomicGroup = NULL,
             quiet = FALSE) {
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
        baseUrl <- pasteURL(
            "ftp.ncbi.nlm.nih.gov", "genomes", "refseq",
            protocol = "ftp"
        )
        if (is.null(taxonomicGroup)) {
            taxonomicGroup <- .matchNcbiTaxonomicGroup(
                organism = organism,
                mode = "refseq"
            )
        }
        url <- pasteURL(
            baseUrl,
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
#' @note Updated 2022-02-08.
#' @noRd
#'
#' @param file `character(1)`.
#' RefSeq GFF file.
#'
#' @return `Seqinfo`.
#'
#' @seealso
#' - `tximeta:::gtf2RefSeq()`.
#'
#' @examples
#' ## RefSeq GRCh38.p12.
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
#'
#' ## RefSeq GRCh38 assembly for pipelines.
#' file <- pasteURL(
#'     "ftp.ncbi.nlm.nih.gov",
#'     "genomes",
#'     "all",
#'     "GCA",
#'     "000",
#'     "001",
#'     "405",
#'     "GCA_000001405.15_GRCh38",
#'     "seqs_for_alignment_pipelines.ucsc_ids",
#'     "GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gff.gz",
#'     protocol = "ftp"
#' )
#' seqinfo <- .getRefSeqSeqinfo(file)
#' print(seqinfo)
.getRefSeqSeqinfo <- function(file) {
    if (
        identical(
            x = "seqs_for_alignment_pipelines.ucsc_ids",
            y = basename(dirname(file))
        ) ||
            grepl(
                pattern = "_full_analysis_set.refseq_annotation",
                x = basename(file),
                fixed = TRUE
            )
    ) {
        ucsc <- TRUE
    } else {
        ucsc <- FALSE
    }
    ## Locate the "*_assembly_report.txt" file from the GFF file path.
    file <- .locateRefSeqAssemblyReport(file)
    file <- .cacheIt(file)
    pattern <- paste0(
        "^([a-z0-9]+_)?",
        "(GC[AF]_[0-9]+\\.[0-9]+)",
        "_(.+)",
        "_assembly_report",
        "(.+ucsc_names)?",
        "\\.txt$"
    )
    assert(isMatchingRegex(x = basename(file), pattern = pattern))
    alert(sprintf(
        "Getting {.cls %s} from {.file %s}.",
        "Seqinfo", basename(file)
    ))
    ## e.g. "GRCh38.p13", which is the format Seqinfo expects.
    ## Refer to GenomeInfoDb documentation for details on NCBI.
    genomeBuild <- sub(
        pattern = pattern,
        replacement = "\\3",
        x = basename(file)
    )
    ## Need to parse the comments in the assembly file to extract the column
    ## names for the data frame. Note that the number of columns differs in
    ## the report file for the "seqs_for_alignment_pipelines.ucsc_ids" assembly.
    lines <- import(con = file, format = "lines", quiet = TRUE)
    comments <- grep(pattern = "^#", x = lines, value = TRUE)
    colnames <- comments[length(comments)]
    assert(isMatchingFixed(x = colnames, pattern = "\t"))
    colnames <- sub(pattern = "^# ", replacement = "", x = colnames)
    colnames <- strsplit(colnames, split = "\t")[[1L]]
    colnames <- camelCase(colnames, strict = TRUE)
    seqnames <- ifelse(
        test = ucsc,
        yes = "ucscStyleName",
        no = "refSeqAccn"
    )
    whatCols <- c(
        "seqnames" = seqnames,
        "seqlengths" = "sequenceLength"
    )
    assert(isSubset(whatCols, colnames))
    ## NOTE `data.table::fread` doesn't currently support comment exclusion,
    ## so we're intentionally using the base R engine for import here instead.
    df <- import(
        con = file,
        format = "tsv",
        colnames = colnames,
        comment = "#",
        engine = "base",
        quiet = FALSE
    )
    df <- df[, whatCols, drop = FALSE]
    df <- df[complete.cases(df), , drop = FALSE]
    seqnames <- df[[whatCols[["seqnames"]]]]
    seqlengths <- df[[whatCols[["seqlengths"]]]]
    assert(
        !anyNA(seqnames),
        !anyNA(seqlengths),
        hasNoDuplicates(seqnames)
    )
    seq <- Seqinfo(
        seqnames = seqnames,
        seqlengths = seqlengths,
        isCircular = NA,
        genome = genomeBuild
    )
    assert(
        is(seq, "Seqinfo"),
        validObject(seq)
    )
    seq
}



#' Locate RefSeq assembly report, from GFF file
#'
#' @note Updated 2022-02-08.
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
#' ## RefSeq assembly for alignment pipelines.
#' file <- pasteURL(
#'     "ftp.ncbi.nlm.nih.gov",
#'     "genomes",
#'     "all",
#'     "GCA",
#'     "000",
#'     "001",
#'     "405",
#'     "GCA_000001405.15_GRCh38",
#'     "seqs_for_alignment_pipelines.ucsc_ids",
#'     "GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gff.gz",
#'     protocol = "ftp"
#' )
#' x <- .locateRefSeqAssemblyReport(file)
#' print(x)
#'
#' ## `downloadRefSeqGenome` convention.
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
        x = basename(file),
        pattern = .gffPatterns[["refseq"]]
    ))
    reportBasename <- sub(
        pattern = paste0(
            "_(genomic|full_analysis_set.refseq_annotation)",
            "\\.(gff|gtf)(\\.gz)?$"
        ),
        replacement = "_assembly_report\\.txt",
        x = basename(file)
    )
    if (isAURL(file)) {
        if (identical(
            x = "seqs_for_alignment_pipelines.ucsc_ids",
            y = basename(dirname(file))
        )) {
            x <- pasteURL(parentDir(file, n = 2L), reportBasename)
        } else {
            x <- pasteURL(dirname(file), reportBasename)
        }
        return(x)
    }
    ## Full local FTP pulldown.
    x <- file.path(dirname(file), reportBasename)
    if (isAFile(x)) {
        return(x)
    }
    ## `downloadRefSeqGenome()` download.
    x <- file.path(parentDir(file, n = 2L), "metadata", reportBasename)
    if (isAFile(x)) {
        return(x)
    }
    abort(sprintf(
        "Failed to locate RefSeq assembly report from {.file %s}.", file
    ))
}
