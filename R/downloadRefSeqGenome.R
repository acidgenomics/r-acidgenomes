#' Download RefSeq reference genome
#'
#' @export
#' @note Updated 2021-01-08.
#'
#' @inheritParams currentGenomeBuild
#' @inheritParams downloadEnsemblGenome
#'
#' @seealso
#' - [RefSeq Genomes FTP server](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/)
#' - [Human Genome Resources at NCBI](https://www.ncbi.nlm.nih.gov/projects/genome/guide/human/)
#' - [Genomes Download (FTP) FAQ](https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/)
#'
#' @examples
#' ## This example is bandwidth intensive.
#' ## > downloadRefSeqGenome(
#' ## >     organism = "Homo sapiens",
#' ## >     taxonomicGroup = "vertebrate_mammalian"
#' ## >     genomeBuild = "GRCh38",
#' ## >     type = "transcriptome",
#' ## >     annotation = "gtf"
#' ## > )
downloadRefSeqGenome <-
    function(
        organism,
        taxonomicGroup = NULL,
        type = c("all", "transcriptome", "genome", "none"),
        annotation = c("all", "gtf", "gff", "none"),
        outputDir = "."
    ) {
        assert(
            isOrganism(organism),
            isString(taxonomicGroup, nullOK = TRUE),
            isString(genomeBuild, nullOK = TRUE),
            isString(outputDir)
        )
        baseURL <- .getRefSeqGenomeURL(
            organism = organism,
            taxonomicGroup = taxonomicGroup,
            quiet = FALSE
        )
        taxonomicGroup <- basename(dirname(baseURL))
        release <- currentRefSeqVersion()
        summary <- .getRefSeqAssemblySummary(baseURL)
        assert(isSubset("ftp_path", names(summary)))
        releaseURL <- summary[["ftp_path"]]
        genomeBuild <- basename(releaseURL)
        assert(isAURL(releaseURL))
        outputDir <- initDir(outputDir)
        type <- match.arg(type)
        annotation <- match.arg(annotation)
        if (type == "none" && annotation == "none") {
            stop("'type' or 'annotation' argument is required.")
        }
        dlList <- list(
            "type" = c(
                "genome" = FALSE,
                "transcriptome" = FALSE
            ),
            "annotation" = c(
                "gff" = FALSE,
                "gtf" = FALSE
            )
        )
        switch(
            EXPR = type,
            "all" = {
                dlList[["type"]][["genome"]] <- TRUE
                dlList[["type"]][["transcriptome"]] <- TRUE
            },
            "genome" = {
                dlList[["type"]][["genome"]] <- TRUE
            },
            "transcriptome" = {
                dlList[["type"]][["transcriptome"]] <- TRUE
            }
        )
        switch(
            EXPR = annotation,
            "all" = {
                dlList[["annotation"]][["gff"]] <- TRUE
                dlList[["annotation"]][["gtf"]] <- TRUE
            },
            "gff" = {
                dlList[["annotation"]][["gff"]] <- TRUE
            },
            "gtf" = {
                dlList[["annotation"]][["gtf"]] <- TRUE
            }
        )
        outputBasename <- kebabCase(tolower(paste(
            organism, genomeBuild, "refseq", release
        )))
        outputDir <- file.path(outputDir, outputBasename)
        h1(sprintf(
            paste(
                "Downloading RefSeq genome for {.emph %s}",
                " %s %d from {.url %s} to {.path %s}."
            ),
            organism, genomeBuild, release,
            releaseURL, outputDir
        ))
        assert(!isADir(outputDir))
        outputDir <- initDir(outputDir)
        out <- list("summary" = summary)
        urls <- c(
            "readme" = pasteURL(releaseURL, "README.txt"),
            "annotationHashes" = pasteURL(releaseURL, "annotation_hashes.txt"),
            "assemblyStatus" = pasteURL(releaseURL, "assembly_status.txt"),
            "md5checksums" = pasteURL(releaseURL, "md5checksums.txt"),
            "assemblyRegions" = pasteURL(
                releaseURL,
                paste0(genomeBuild, "_assembly_regions.txt")
            ),
            "assemblyReport" = pasteURL(
                releaseURL,
                paste0(genomeBuild, "_assembly_report.txt")
            ),
            "assemblyStats" = pasteURL(
                releaseURL,
                paste0(genomeBuild, "_assembly_stats.txt")
            )
        )
        out[["metadata"]] <- .downloadURLs(urls = urls, outputDir = outputDir)
        ## Now ready to download individual genome files.
        args <- list(
            genomeBuild = genomeBuild,
            releaseURL = releaseURL,
            outputDir = outputDir
        )
        if (isTRUE(dlList[["type"]][["genome"]])) {
            out[["type"]][["genome"]] <-
                do.call(what = .downloadRefSeqFASTA, args = args)
        }
        if (isTRUE(dlList[["type"]][["transcriptome"]])) {
            out[["type"]][["transcriptome"]] <-
                do.call(what = .downloadRefSeqTranscriptomeFASTA, args = args)
        }
        args <- c(args, release = release)
        if (isTRUE(dlList[["annotation"]][["gtf"]])) {
            out[["annotation"]][["gtf"]] <-
                do.call(what = .downloadRefSeqGTF, args = args)
        }
        if (isTRUE(dlList[["annotation"]][["gff"]])) {
            out[["annotation"]][["gff"]] <-
                do.call(what = .downloadRefSeqGFF, args = args)
        }
        out[["args"]] <- args
        out[["call"]] <- match.call()
        saveRDS(
            object = sessionInfo(),
            file = file.path(outputDir, "sessionInfo.rds")
        )
        alertSuccess(sprintf(
            "RefSeq genome downloaded successfully to {.path %s}.",
            outputDir
        ))
        invisible(out)
    }



## Updated 2021-01-08.
.downloadRefSeqGFF <-
    function(
        genomeBuild,
        releaseURL,
        outputDir
    ) {
        urls <- c(
            "gff" = pasteURL(
                releaseURL,
                paste0(genomeBuild, "_genomic.gff.gz")
            )
        )
        .downloadURLs(urls = urls, outputDir = file.path(outputDir, "gff"))
    }



## Updated 2021-01-08.
.downloadRefSeqGTF <-
    function(
        genomeBuild,
        releaseURL,
        outputDir
    ) {
        urls <- c(
            "gtf" = pasteURL(
                releaseURL,
                paste0(genomeBuild, "_genomic.gtf.gz")
            )
        )
        .downloadURLs(urls = urls, outputDir = file.path(outputDir, "gtf"))
    }



## Updated 2021-01-08.
.downloadRefSeqGenomeFASTA <-
    function(
        genomeBuild,
        releaseURL,
        outputDir
    ) {
        urls <- c(
            "fasta" = pasteURL(
                releaseURL,
                paste0(genomeBuild, "_genomic.fna.gz")
            )
        )
        .downloadURLs(urls = urls, outputDir = file.path(outputDir, "genome"))
    }



## Updated 2021-01-08.
.downloadRefSeqTranscriptomeFASTA <-
    function(
        genomeBuild,
        releaseURL,
        outputDir
    ) {
        out <- list()
        urls <- c(
            "fasta" = pasteURL(
                releaseURL,
                paste0(genomeBuild, "_rna.fna.gz")
            )
        )
        out[["fasta"]] <- .downloadURLs(
            urls = urls,
            outputDir = file.path(outputDir, "transcriptome")
        )
        fastaFile <- out[["fasta"]]
        assert(isAFile(fastaFile))
        ## FIXME THIS ISNT SUPPORTED YET IN THE PACKAGE.
        out[["tx2gene"]] <- makeTx2GeneFileFromFASTA(
            file = fastaFile,
            source = "refseq"
        )
        invisible(out)
    }
