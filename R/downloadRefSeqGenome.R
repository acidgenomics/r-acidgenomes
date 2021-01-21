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
#' @note Updated 2021-01-21.
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
        args <- list(
            "genomeBuild" = genomeBuild,
            "outputDir" = outputDir,
            "releaseURL" = releaseURL
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
        info[["call"]] <- match.call()
        info[["sessionInfo"]] <- sessionInfo()
        saveRDS(object = info, file = file.path(outputDir, "metadata.rds"))
        alertSuccess(sprintf(
            "RefSeq genome downloaded successfully to {.path %s}.",
            outputDir
        ))
        invisible(info)
    }



## Updated 2021-01-20.
.downloadRefSeqAnnotation <-
    function(
        genomeBuild,
        outputDir,
        releaseURL
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
            outputDir = file.path(outputDir, "annotation")
        )
        ## Create GFF and GTF symlinks.
        if (!isWindows()) {
            gffFile <- files[["gff"]]
            assert(isAFile(gffFile))
            gffSymlink <- file.path(outputDir, "annotation.gff3.gz")
            file.symlink(from = gffFile, to = gffSymlink)
            files[["gffSymlink"]] <- gffSymlink
            ## Generate GTF symlink.
            gtfFile <- files[["gtf"]]
            assert(isAFile(gtfFile))
            gtfSymlink <- file.path(
                outputDir,
                paste0("annotation.", fileExt(gtfFile))
            )
            file.symlink(from = gtfFile, to = gtfSymlink)
            files[["gtfSymlink"]] <- gtfSymlink
        }
        invisible(list("files" = files, "urls" = urls))
    }



## Updated 2021-01-20.
.downloadRefSeqGenome <-
    function(
        genomeBuild,
        outputDir,
        releaseURL
    ) {
        urls <- c(
            "fasta" = pasteURL(
                releaseURL,
                paste0(genomeBuild, "_genomic.fna.gz")
            )
        )
        files <- .downloadURLs(
            urls = urls,
            outputDir = file.path(outputDir, "genome")
        )
        ## Create FASTA symlink.
        if (!isWindows()) {
            fastaFile <- files[["fasta"]]
            assert(isAFile(fastaFile))
            fastaSymlink <- file.path(outputDir, "genome.fa.gz")
            file.symlink(from = fastaFile, to = fastaSymlink)
            files[["fastaSymlink"]] <- fastaSymlink
        }
        invisible(list("files" = files, "urls" = urls))
    }



## Updated 2021-01-20.
.downloadRefSeqMetadata <-
    function(
        genomeBuild,
        outputDir,
        releaseURL
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
            outputDir = file.path(outputDir, "metadata")
        )
        invisible(list("files" = files, "urls" = urls))
    }


## Updated 2021-01-20.
.downloadRefSeqTranscriptome <-
    function(
        genomeBuild,
        outputDir,
        releaseURL
    ) {
        urls <- c(
            "fasta" = pasteURL(
                releaseURL,
                paste0(genomeBuild, "_rna.fna.gz")
            )
        )
        files <- .downloadURLs(
            urls = urls,
            outputDir = file.path(outputDir, "transcriptome")
        )
        ## Create FASTA symlink.
        if (!isWindows()) {
            fastaFile <- files[["fasta"]]
            assert(isAFile(fastaFile))
            fastaSymlink <- file.path(outputDir, "transcriptome.fa.gz")
            file.symlink(from = fastaFile, to = fastaSymlink)
            files[["fastaSymlink"]] <- fastaSymlink
        }
        invisible(list("files" = files, "urls" = urls))
    }
