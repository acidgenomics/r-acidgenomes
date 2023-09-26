## nolint start
#' Download UCSC reference genome
#'
#' @export
#' @note Updated 2023-09-26.
#'
#' @section Genome:
#'
#' - `<GENOME_BUILD>.chrom.sizes`: Two-column tab-separated text file containing
#' assembly sequence names and sizes.
#' - `<GENOME_BUILD>.chromAlias.txt`: Sequence name alias file, one line for
#' each sequence name.  First column is sequence name followed by tab
#' separated alias names.
#'
#' @section Transcriptome:
#'
#' - `mrna.fa.gz`: Human mRNA from GenBank. This sequence data is updated
#' regularly via automatic GenBank updates.
#' - `refMrna.fa.gz`: RefSeq mRNA from the same species as the genome.
#' This sequence data is updated regularly via automatic GenBank updates.
#'
#' @section Gene annotations:
#'
#' This directory contains GTF files for the main gene transcript sets where
#' available. They are sourced from the following gene model tables:
#' `knownGene` (GENCODE) and `ncbiRefSeq` (NCBI RefSeq).
#'
#' @inheritParams currentGenomeBuild
#' @inheritParams downloadEnsemblGenome
#'
#' @param genomeBuild `character(1)`.
#' UCSC genome build assembly name (e.g. `"hg38"`).
#' If set `NULL`, defauls to the most recent build available.
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
    function(organism,
             genomeBuild = NULL,
             outputDir = getwd(),
             cache = FALSE) {
        assert(
            isOrganism(organism),
            isString(genomeBuild, nullOK = TRUE),
            isString(outputDir),
            isFlag(cache)
        )
        organism <- match.arg(
            arg = organism,
            choices = c("Homo sapiens", "Mus musculus")
        )
        outputDir <- initDir(outputDir)
        ## Homo sapiens is now defaulting to "hs1", which is the experimental
        ## new T2T-CHM13 assembly. We're not quite ready to support this yet, so
        ## keeping pinned to hg38.
        if (is.null(genomeBuild)) {
            genomeBuild <- switch(
                EXPR = organism,
                "Homo sapiens" = "hg38",
                currentUCSCGenomeBuild(organism)
            )
        }
        ## UCSC is updated on a rolling schedule, so use today's date as a
        ## substitution for release number.
        release <- Sys.Date()
        ## UCSC requests downloads over FTP instead of HTTPS when possible,
        ## so respecting their wishes here.
        baseUrl <- pasteURL(
            "hgdownload.soe.ucsc.edu",
            "goldenPath",
            protocol = "ftp"
        )
        releaseUrl <- pasteURL(baseUrl, genomeBuild, "bigZips")
        outputBasename <- kebabCase(tolower(paste(
            organism, genomeBuild, "ucsc", release
        )))
        outputDir <- file.path(outputDir, outputBasename)
        h1(sprintf(
            paste(
                "Downloading UCSC genome for {.emph %s}",
                "%s from {.url %s} to {.path %s}."
            ),
            organism, genomeBuild,
            releaseUrl, outputDir
        ))
        assert(
            !isADir(outputDir),
            msg = sprintf("Genome exists at '%s'.", outputDir)
        )
        outputDir <- initDir(outputDir)
        args <- list(
            "outputDir" = outputDir,
            "releaseUrl" = releaseUrl,
            "cache" = cache
        )
        info <- list()
        info[["date"]] <- Sys.Date()
        info[["metadata"]] <-
            do.call(what = .downloadUcscMetadata, args = args)
        info[["transcriptome"]] <-
            do.call(what = .downloadUcscTranscriptome, args = args)
        args <- append(x = args, values = list("genomeBuild" = genomeBuild))
        info[["genome"]] <-
            do.call(what = .downloadUcscGenome, args = args)
        info[["annotation"]] <-
            do.call(what = .downloadUcscAnnotation, args = args)
        info[["args"]] <- args
        info[["call"]] <- tryCatch(
            expr = standardizeCall(),
            error = function(e) {
                NULL
            }
        )
        info[["sessionInfo"]] <- sessionInfo()
        saveRDS(object = info, file = file.path(outputDir, "metadata.rds"))
        alertSuccess(sprintf(
            "UCSC genome downloaded successfully to {.path %s}.",
            outputDir
        ))
        invisible(info)
    }



## Updated 2023-09-26.
.downloadUcscAnnotation <-
    function(genomeBuild,
             outputDir,
             releaseUrl,
             cache) {
        genesUrl <- pasteURL(releaseUrl, "genes")
        urls <- c(
            "readme" = pasteURL(genesUrl, "README.txt"),
            "knownGene" = pasteURL(
                genesUrl,
                paste0(genomeBuild, ".knownGene.gtf.gz")
            ),
            "ncbiRefSeq" = pasteURL(
                genesUrl,
                paste0(genomeBuild, ".ncbiRefSeq.gtf.gz")
            )
        )
        files <- .downloadUrls(
            urls = urls,
            outputDir = file.path(outputDir, "annotation"),
            cache = cache
        )
        gtfFile <- files[["knownGene"]]
        ## Create relative path symlink.
        if (!isWindows() && requireNamespace("withr", quietly = TRUE)) {
            gtfRelativeFile <- sub(
                pattern = paste0("^", outputDir, "/"),
                replacement = "",
                x = gtfFile
            )
            gtfSymlink <- paste0("annotation.", fileExt(gtfFile))
            withr::with_dir(
                new = outputDir,
                code = {
                    file.symlink(from = gtfRelativeFile, to = gtfSymlink)
                }
            )
            files[["gtfSymlink"]] <- gtfSymlink
        }
        ## Save genomic ranges.
        genes <- makeGRangesFromGFF(gtfFile, level = "genes")
        transcripts <- makeGRangesFromGFF(gtfFile, level = "transcripts")
        saveRDS(
            object = genes,
            file = file.path(outputDir, "genes.rds")
        )
        saveRDS(
            object = transcripts,
            file = file.path(outputDir, "transcripts.rds")
        )
        ## Save transcript-to-gene mappings.
        t2g <- makeTxToGeneFromGFF(gtfFile)
        saveRDS(object = t2g, file = file.path(outputDir, "tx2gene.rds"))
        t2gFile <- export(
            object = t2g,
            con = file.path(outputDir, "tx2gene.csv.gz")
        )
        files[["tx2gene"]] <- t2gFile
        invisible(list("files" = files, "urls" = urls))
    }



## Note that both hg38 and hg19 support "latest/" subdirectory.
## Updated 2022-05-24.
.downloadUcscGenome <-
    function(genomeBuild,
             outputDir,
             releaseUrl,
             cache) {
        isHuman <- grepl(pattern = "^hg[0-9]+$", x = genomeBuild)
        latestUrl <- ifelse(
            test = isHuman,
            yes = pasteURL(releaseUrl, "latest"),
            no = releaseUrl
        )
        urls <- c(
            "chromAlias" = pasteURL(
                releaseUrl,
                paste0(genomeBuild, ".chromAlias.txt")
            ),
            "chromSizes" = pasteURL(
                latestUrl,
                paste0(genomeBuild, ".chrom.sizes")
            ),
            "fasta" = pasteURL(
                latestUrl,
                paste0(genomeBuild, ".fa.gz")
            ),
            "md5sum" = pasteURL(latestUrl, "md5sum.txt")
        )
        if (identical(genomeBuild, "hg38")) {
            urls <- c(
                urls,
                "latestVersion" = pasteURL(
                    releaseUrl, "latest", "LATEST_VERSION"
                )
            )
        }
        files <- .downloadUrls(
            urls = urls,
            outputDir = file.path(outputDir, "genome"),
            cache = cache
        )
        ## Create relative path symlink.
        if (!isWindows() && requireNamespace("withr", quietly = TRUE)) {
            fastaFile <- files[["fasta"]]
            fastaRelativeFile <- sub(
                pattern = paste0("^", outputDir, "/"),
                replacement = "",
                x = fastaFile
            )
            fastaSymlink <- paste0("genome.", fileExt(fastaFile))
            withr::with_dir(
                new = outputDir,
                code = {
                    file.symlink(from = fastaRelativeFile, to = fastaSymlink)
                }
            )
            files[["fastaSymlink"]] <- fastaSymlink
        }
        invisible(list("files" = files, "urls" = urls))
    }



## Updated 2021-08-03.
.downloadUcscMetadata <-
    function(outputDir,
             releaseUrl,
             cache) {
        urls <- c("readme" = pasteURL(releaseUrl, "README.txt"))
        files <- .downloadUrls(
            urls = urls,
            outputDir = file.path(outputDir, "metadata"),
            cache = cache
        )
        invisible(list("files" = files, "urls" = urls))
    }



## Updated 2022-05-24.
.downloadUcscTranscriptome <-
    function(outputDir,
             releaseUrl,
             cache) {
        urls <- c(
            "mrna" = pasteURL(releaseUrl, "mrna.fa.gz"),
            "mrnaChecksum" = pasteURL(releaseUrl, "mrna.fa.gz.md5"),
            "refMrna" = pasteURL(releaseUrl, "refMrna.fa.gz"),
            "refMrnaChecksum" = pasteURL(releaseUrl, "refMrna.fa.gz.md5")
        )
        files <- .downloadUrls(
            urls = urls,
            outputDir = file.path(outputDir, "metadata"),
            cache = cache
        )
        ## Create relative path symlink.
        if (!isWindows() && requireNamespace("withr", quietly = TRUE)) {
            fastaFile <- files[["mrna"]]
            fastaRelativeFile <- sub(
                pattern = paste0("^", outputDir, "/"),
                replacement = "",
                x = fastaFile
            )
            fastaSymlink <- paste0("transcriptome.", fileExt(fastaFile))
            withr::with_dir(
                new = outputDir,
                code = {
                    file.symlink(from = fastaRelativeFile, to = fastaSymlink)
                }
            )
            files[["fastaSymlink"]] <- fastaSymlink
        }
        invisible(list("files" = files, "urls" = urls))
    }
