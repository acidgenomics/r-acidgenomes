## FIXME Need to rework the `genomeBuild` documentation here.

## nolint start

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
#' @note Updated 2021-08-03.
#'
#' @inheritParams currentGenomeBuild
#' @inheritParams downloadEnsemblGenome
#'
#' @return Invisible `list`.
#'
#' @seealso
#' - [Human Genome Resources at NCBI](https://www.ncbi.nlm.nih.gov/projects/genome/guide/human/)
#' - [RefSeq Genomes FTP server](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/)
#' - [Genomes Download (FTP) FAQ](https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/)
#' - [bcbio hg38 reference genome](https://steinbaugh.com/posts/bcbio-hg38.html)
#' - [Heng Li: Which human reference genome to use?](http://lh3.github.io/2017/11/13/which-human-reference-genome-to-use)
#' - [GRCh38 assembly for alignment pipelines](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/)
#' - [UCSC hg38 bigZips](https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/)
#'
#' @examples
#' ## This example is bandwidth intensive.
#' ## > downloadRefSeqGenome(
#' ## >     organism = "Homo sapiens",
#' ## >     taxonomicGroup = "vertebrate_mammalian",
#' ## >     genomeBuild = "GCF_000001405.39_GRCh38.p12",
#' ## > )

## nolint end

downloadRefSeqGenome <-
    function(
        organism,
        taxonomicGroup = NULL,
        genomeBuild = NULL,
        outputDir = ".",
        cache = TRUE
    ) {
        assert(
            isOrganism(organism),
            isString(taxonomicGroup, nullOK = TRUE),
            isString(genomeBuild, nullOK = TRUE),
            isString(outputDir),
            isFlag(cache)
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
        assert(
            !isADir(outputDir),
            msg = sprintf("Genome exists at '%s'.", outputDir)
        )
        outputDir <- initDir(outputDir)
        args <- list(
            "genomeBuild" = genomeBuild,
            "outputDir" = outputDir,
            "releaseURL" = releaseURL,
            "cache" = cache
        )
        info <- list()
        info[["date"]] <- Sys.Date()
        info[["metadata"]] <-
            do.call(what = .downloadRefSeqMetadata, args = args)
        info[["genome"]] <-
            do.call(what = .downloadRefSeqGenome, args = args)
        info[["transcriptome"]] <-
            do.call(what = .downloadRefSeqTranscriptome, args = args)
        info[["annotation"]] <-
            do.call(what = .downloadRefSeqAnnotation, args = args)
        info[["args"]] <- args
        info[["call"]] <- tryCatch(
            expr = standardizeCall(),
            error = function(e) NULL
        )
        info[["sessionInfo"]] <- sessionInfo()
        saveRDS(object = info, file = file.path(outputDir, "metadata.rds"))
        alertSuccess(sprintf(
            "RefSeq genome downloaded successfully to {.path %s}.",
            outputDir
        ))
        invisible(info)
    }



## Updated 2021-08-03.
.downloadRefSeqAnnotation <-
    function(
        genomeBuild,
        outputDir,
        releaseURL,
        cache
    ) {
        urls <- c(
            "gff" = pasteURL(
                releaseURL,
                paste0(genomeBuild, "_genomic.gff.gz")
            ),
            "gtf" = pasteURL(
                releaseURL,
                paste0(genomeBuild, "_genomic.gtf.gz")
            )
        )
        files <- .downloadURLs(
            urls = urls,
            outputDir = file.path(outputDir, "annotation"),
            cache = cache
        )
        gffFile <- files[["gff"]]
        gtfFile <- files[["gtf"]]
        ## Create GFF and GTF symlinks.
        if (!isWindows()) {
            wd <- getwd()
            setwd(outputDir)
            gffRelativeFile <- sub(
                pattern = paste0("^", outputDir, "/"),
                replacement = "",
                x = gffFile
            )
            gtfRelativeFile <- sub(
                pattern = paste0("^", outputDir, "/"),
                replacement = "",
                x = gtfFile
            )
            assert(
                isAFile(gffFile),
                isAFile(gffRelativeFile),
                isAFile(gtfFile),
                isAFile(gtfRelativeFile)
            )
            gffSymlink <- "annotation.gff3.gz"
            gtfSymlink <- "annotation.gtf.gz"
            file.symlink(from = gffRelativeFile, to = gffSymlink)
            file.symlink(from = gtfRelativeFile, to = gtfSymlink)
            files[["gffSymlink"]] <- gffSymlink
            files[["gtfSymlink"]] <- gtfSymlink
            setwd(wd)
        }
        ## Save genomic ranges.
        genes <- makeGRangesFromGFF(gffFile, level = "genes")
        transcripts <- makeGRangesFromGFF(gffFile, level = "transcripts")
        saveRDS(
            object = genes,
            file = file.path(outputDir, "genes.rds")
        )
        saveRDS(
            object = transcripts,
            file = file.path(outputDir, "transcripts.rds")
        )
        ## Save transcript-to-gene mappings.
        tx2gene <- makeTx2GeneFromGFF(file = gffFile)
        saveRDS(object = tx2gene, file = file.path(outputDir, "tx2gene.rds"))
        tx2geneFile <- export(
            object = tx2gene,
            file = file.path(outputDir, "tx2gene.csv.gz")
        )
        files[["tx2gene"]] <- tx2geneFile
        invisible(list("files" = files, "urls" = urls))
    }



## Updated 2021-08-03.
.downloadRefSeqGenome <-
    function(
        genomeBuild,
        outputDir,
        releaseURL,
        cache
    ) {
        urls <- c(
            "fasta" = pasteURL(
                releaseURL,
                paste0(genomeBuild, "_genomic.fna.gz")
            )
        )
        files <- .downloadURLs(
            urls = urls,
            outputDir = file.path(outputDir, "genome"),
            cache = cache
        )
        ## Create FASTA symlink.
        if (!isWindows()) {
            wd <- getwd()
            setwd(outputDir)
            fastaFile <- files[["fasta"]]
            fastaRelativeFile <- sub(
                pattern = paste0("^", outputDir, "/"),
                replacement = "",
                x = fastaFile
            )
            assert(
                isAFile(fastaFile),
                isAFile(fastaRelativeFile)
            )
            fastaSymlink <- "genome.fa.gz"
            file.symlink(from = fastaRelativeFile, to = fastaSymlink)
            files[["fastaSymlink"]] <- fastaSymlink
            setwd(wd)
        }
        invisible(list("files" = files, "urls" = urls))
    }



## Updated 2021-08-03.
.downloadRefSeqMetadata <-
    function(
        genomeBuild,
        outputDir,
        releaseURL,
        cache
    ) {
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
        files <- .downloadURLs(
            urls = urls,
            outputDir = file.path(outputDir, "metadata"),
            cache = cache
        )
        invisible(list("files" = files, "urls" = urls))
    }



## Updated 2021-08-03.
.downloadRefSeqTranscriptome <-
    function(
        genomeBuild,
        outputDir,
        releaseURL,
        cache
    ) {
        urls <- c(
            "fasta" = pasteURL(
                releaseURL,
                paste0(genomeBuild, "_rna.fna.gz")
            )
        )
        files <- .downloadURLs(
            urls = urls,
            outputDir = file.path(outputDir, "transcriptome"),
            cache = cache
        )
        ## Create FASTA symlink.
        if (!isWindows()) {
            wd <- getwd()
            setwd(outputDir)
            fastaFile <- files[["fasta"]]
            fastaRelativeFile <- sub(
                pattern = paste0("^", outputDir, "/"),
                replacement = "",
                x = fastaFile
            )
            assert(
                isAFile(fastaFile),
                isAFile(fastaRelativeFile)
            )
            fastaSymlink <- "transcriptome.fa.gz"
            file.symlink(from = fastaRelativeFile, to = fastaSymlink)
            files[["fastaSymlink"]] <- fastaSymlink
            setwd(wd)

        }
        invisible(list("files" = files, "urls" = urls))
    }
