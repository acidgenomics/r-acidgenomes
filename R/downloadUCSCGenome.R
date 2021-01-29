## nolint start

#' Download UCSC reference genome
#'
#' @export
#' @note Updated 2021-01-29.
#'
#' @section Genome:
#'
#' - `<GENOME_BUILD>.chrom.sizes`: Two-column tab-separated text file containing
#'   assembly sequence names and sizes.
#'  - `<GENOME_BUILD>.chromAlias.txt`: Sequence name alias file, one line for
#'    each sequence name.  First column is sequence name followed by tab
#'    separated alias names.
#'
#' @section Transcriptome:
#'
#' - `mrna.fa.gz`: Human mRNA from GenBank. This sequence data is updated
#'   regularly via automatic GenBank updates.
#' - `refMrna.fa.gz`: RefSeq mRNA from the same species as the genome.
#'   This sequence data is updated regularly via automatic GenBank updates.
#'
#' @section Gene annotations:
#' :
#' This directory contains GTF files for the main gene transcript sets where
#' available. They are sourced from the following gene model tables:
#' ncbiRefSeq, refGene, ensGene, knownGene.
#'
#' @inheritParams currentGenomeBuild
#' @inheritParams downloadEnsemblGenome
#'
#' @return Invisible `list`.
#'
#' @seealso
#' - [hg38 pinned analysis set (for NGS pipelines)](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/).
#'
#' @examples
#' ## This example is bandwidth intensive.
#' ## > downloadUCSCGenome(organism = "Homo sapiens")

## nolint end

downloadUCSCGenome <-
    function(
        organism,
        genomeBuild = NULL,
        outputDir = "."
    ) {
        assert(
            isOrganism(organism),
            isString(genomeBuild, nullOK = TRUE),
            isString(outputDir)
        )
        outputDir <- initDir(outputDir)
        if (is.null(genomeBuild)) {
            genomeBuild <- currentUCSCGenomeBuild(organism)
        }
        ## UCSC is updated on a rolling schedule, so use today's date as a
        ## substitution for release number.
        release <- Sys.Date()
        ## Consider downloading over rsync (preferred by UCSC) instead of FTP
        ## in a future update. This requires updating internal downloader to
        ## support rsync, which is currently only useful here.
        ## > protocol <- ifelse(
        ## >     test = isTRUE(isSystemCommand("rsync")),
        ## >     yes = "rsync",
        ## >     no = "ftp"
        ## > )
        baseURL <- pasteURL(
            "hgdownload.soe.ucsc.edu",
            "goldenPath",
            protocol = "ftp"
        )
        releaseURL <- pasteURL(baseURL, genomeBuild, "bigZips")
        outputBasename <- kebabCase(tolower(paste(
            organism, genomeBuild, "ucsc", release
        )))
        outputDir <- file.path(outputDir, outputBasename)
        h1(sprintf(
            paste(
                "Downloading Ensembl genome for {.emph %s}",
                " %s from {.url %s} to {.path %s}."
            ),
            organism, genomeBuild,
            releaseURL, outputDir
        ))
        assert(!isADir(outputDir))
        outputDir <- initDir(outputDir)
        args <- list(
            "outputDir" = outputDir,
            "releaseURL" = releaseURL
        )
        info <- list()
        info[["date"]] <- Sys.Date()
        info[["metadata"]] <-
            do.call(what = .downloadUCSCMetadata, args = args)
        info[["transcriptome"]] <-
            do.call(what = .downloadUCSCTranscriptome, args = args)
        args <- append(x = args, values = list("genomeBuild" = genomeBuild))
        info[["genome"]] <-
            do.call(what = .downloadUCSCGenome, args = args)
        info[["annotation"]] <-
            do.call(what = .downloadUCSCAnnotation, args = args)
        info[["args"]] <- args
        info[["call"]] <- match.call()
        info[["sessionInfo"]] <- sessionInfo()
        saveRDS(object = info, file = file.path(outputDir, "metadata.rds"))
        alertSuccess(sprintf(
            "UCSC genome downloaded successfully to {.path %s}.",
            outputDir
        ))
        invisible(info)
    }



## Updated 2021-01-29.
.downloadUCSCAnnotation <-
    function(
        genomeBuild,
        outputDir,
        releaseURL
    ) {
        genesURL <- pasteURL(releaseURL, "genes")
        urls <- c(
            "readme" = pasteURL(genesURL, "README.txt"),
            "ensGene" = pasteURL(
                genesURL,
                paste0(genomeBuild, ".ensGene.gtf.gz")
            ),
            "knownGene" = pasteURL(
                genesURL,
                paste0(genomeBuild, ".knownGene.gtf.gz")
            ),
            "ncbiRefSeq" = pasteURL(
                genesURL,
                paste0(genomeBuild, ".ncbiRefSeq.gtf.gz")
            ),
            "refGene" = pasteURL(
                genesURL,
                paste0(genomeBuild, ".refGene.gtf.gz")
            )
        )
        files <- .downloadURLs(
            urls = urls,
            outputDir = file.path(outputDir, "annotation")
        )
        ## Create symlink.
        gtfFile <- files[["ensGene"]]
        if (!isWindows()) {
            assert(isAFile(gtfFile))
            gtfSymlink <- file.path(
                outputDir,
                paste0("annotation.", fileExt(gtfFile))
            )
            file.symlink(from = gtfFile, to = gtfSymlink)
            files[["gtfSymlink"]] <- gtfSymlink
        }
        tx2geneFile <- makeTx2GeneFileFromGFF(file = gtfFile)
        files[["tx2gene"]] <- tx2geneFile
        invisible(list("files" = files, "urls" = urls))
    }



## Note that both hg38 and hg19 support "latest/" subdirectory.
## Updated 2021-01-21.
.downloadUCSCGenome <-
    function(
        genomeBuild,
        outputDir,
        releaseURL
    ) {
        isHuman <- grepl(pattern = "^hg[0-9]+$", x = genomeBuild)
        latestURL <- ifelse(
            test = isHuman,
            yes = pasteURL(releaseURL, "latest"),
            no = releaseURL
        )
        urls <- c(
            "chromAlias" = pasteURL(
                releaseURL,
                paste0(genomeBuild, ".chromAlias.txt")
            ),
            "chromSizes" = pasteURL(
                latestURL,
                paste0(genomeBuild, ".chrom.sizes")
            ),
            "fasta" = pasteURL(
                latestURL,
                paste0(genomeBuild, ".fa.gz")
            ),
            "md5sum" = pasteURL(latestURL, "md5sum.txt")
        )
        if (identical(genomeBuild, "hg38")) {
            urls <- c(
                urls,
                "latestVersion" = pasteURL(
                    releaseURL, "latest", "LATEST_VERSION"
                )
            )
        }
        files <- .downloadURLs(
            urls = urls,
            outputDir = file.path(outputDir, "genome")
        )
        ## Create FASTA symlink.
        if (!isWindows()) {
            fastaFile <- files[["fasta"]]
            assert(isAFile(fastaFile))
            fastaSymlink <- file.path(
                outputDir,
                paste0("genome.", fileExt(fastaFile))
            )
            file.symlink(from = fastaFile, to = fastaSymlink)
            files[["fastaSymlink"]] <- fastaSymlink
        }
        invisible(list("files" = files, "urls" = urls))
    }



## Updated 2021-01-21.
.downloadUCSCMetadata <-
    function(
        outputDir,
        releaseURL
    ) {
        urls <- c("readme" = pasteURL(releaseURL, "README.txt"))
        files <- .downloadURLs(
            urls = urls,
            outputDir = file.path(outputDir, "metadata")
        )
        invisible(list("files" = files, "urls" = urls))
    }



## Updated 2021-01-21.
.downloadUCSCTranscriptome <-
    function(
        outputDir,
        releaseURL
    ) {
        urls <- c(
            "mrna" = pasteURL(releaseURL, "mrna.fa.gz"),
            "mrnaChecksum" = pasteURL(releaseURL, "mrna.fa.gz.md5"),
            "refMrna" = pasteURL(releaseURL, "refMrna.fa.gz"),
            "refMrnaChecksum" = pasteURL(releaseURL, "refMrna.fa.gz.md5")
        )
        files <- .downloadURLs(
            urls = urls,
            outputDir = file.path(outputDir, "metadata")
        )
        ## Create symlink.
        if (!isWindows()) {
            fastaFile <- files[["mrna"]]
            fastaSymlink <- file.path(
                outputDir,
                paste0("transcriptome.", fileExt(fastaFile))
            )
            file.symlink(from = fastaFile, to = fastaSymlink)
            files[["fastaSymlink"]] <- fastaSymlink
        }
        invisible(list("files" = files, "urls" = urls))
    }
