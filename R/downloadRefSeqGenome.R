## FIXME NEED TO GET TRANSCRIPTS TO FLAT GRANGES, THEN WE CAN ADD TX2GENE SUPPORT?
## Alternatively, don't output a tx2gene for this automatically?



#' Download RefSeq reference genome
#'
#' @section Stable release:
#' The latest assembly defined under the "release/" subdirectory is not
#' considered "stable" by the RefSeq team. It is considered good practice to use
#' a genome build one version behind as a stable release
#' (e.g. "GCF_000001405.38_GRCh38.p12" instead of current
#' "GCF_000001405.39_GRCh38.p13" build).
#'
#' @export
#' @note Updated 2021-01-14.
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
#' ## >     taxonomicGroup = "vertebrate_mammalian",
#' ## >     genomeBuild = "GCF_000001405.39_GRCh38.p12",
#' ## >     type = "transcriptome",
#' ## >     annotation = "gtf"
#' ## > )
downloadRefSeqGenome <-
    function(
        organism,
        taxonomicGroup = NULL,
        genomeBuild = NULL,
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
        release <- currentRefSeqVersion()
        baseURL <- .getRefSeqGenomeURL(
            organism = organism,
            taxonomicGroup = taxonomicGroup,
            quiet = FALSE
        )
        if (is.null(taxonomicGroup)) {
            taxonomicGroup <- basename(dirname(baseURL))
        }
        if (is.null(genomeBuild)) {
            genomeBuild <- currentRefSeqGenomeBuild(
                organism = organism,
                taxonomicGroup = taxonomicGroup
            )
        }
        releaseURL <- pasteURL(baseURL, "all_assembly_versions", genomeBuild)
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
        ## Export transcript-to-gene mappings.
        if (isAFile(out[["annotation"]][["gtf"]])) {
            makeTx2GeneFileFromGTF(file = out[["annotation"]][["gtf"]])
        } else if (isAFile(out[["annotation"]][["gff"]])) {
            makeTx2GeneFileFromGFF(file = out[["annotation"]][["gff"]])
        }

        ## FIXME NEED TO RETHINK TX2GENE HANDLING HERE...

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
        urls <- c(
            "fasta" = pasteURL(
                releaseURL,
                paste0(genomeBuild, "_rna.fna.gz")
            )
        )
        .downloadURLs(
            urls = urls,
            outputDir = file.path(outputDir, "transcriptome")
        )
        invisible(out)
    }
