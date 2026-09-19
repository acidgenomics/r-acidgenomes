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
#' @note Updated 2023-04-14.
#'
#' @inheritParams currentGenomeBuild
#' @inheritParams downloadEnsemblGenome
#'
#' @param genomeBuild `character(1)`.
#' RefSeq genome build assembly name (e.g. `"GCF_000001405.39_GRCh38.p12"`).
#' If set `NULL`, defauls to the most recent build available.
#'
#' @param analysisSet `logical(1)`.
#' Download the NCBI "no-alt analysis set" FASTA alongside the standard
#' genome. This ALT-contig-free FASTA is recommended for aligners such as
#' STAR that do not natively handle ALT contigs. Currently supported for
#' GRCh38 and GRCm39. Defaults to `FALSE`.
#'
#' @return Invisible `list`.
#'
#' @seealso
#' - [Human Genome Resources at NCBI](https://www.ncbi.nlm.nih.gov/projects/genome/guide/human/)
#' - [RefSeq Genomes FTP server](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/)
#' - [Genomes Download (FTP) FAQ](https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/)
#' - [GenBank annotating genomes with GFF3 or GTF files](https://www.ncbi.nlm.nih.gov/genbank/genomes_gff/)
#' - [db_xrefs](https://www.ncbi.nlm.nih.gov/genbank/collab/db_xref/)
#' - [bcbio hg38 reference genome](https://steinbaugh.com/posts/bcbio-hg38.html)
#' - [Heng Li: Which human reference genome to use?](https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use)
#' - [GRCh38 assembly for alignment pipelines](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/)
#' - [UCSC hg38 bigZips](https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/)
#'
#' @examples
#' ## This example is bandwidth intensive.
#' ## > downloadRefseqGenome(
#' ## >     organism = "Homo sapiens",
#' ## >     taxonomicGroup = "vertebrate_mammalian",
#' ## >     genomeBuild = "GCF_000001405.39_GRCh38.p12",
#' ## > )
## nolint end
downloadRefseqGenome <-
    function(
        organism,
        taxonomicGroup = NULL,
        genomeBuild = NULL,
        outputDir = getwd(),
        cache = FALSE,
        analysisSet = FALSE
    ) {
        assert(
            isOrganism(organism),
            isString(taxonomicGroup, nullOk = TRUE),
            isString(genomeBuild, nullOk = TRUE),
            isString(outputDir),
            isFlag(cache),
            isFlag(analysisSet)
        )
        release <- currentRefseqVersion()
        baseUrl <- .getRefSeqGenomeUrl(
            organism = organism,
            taxonomicGroup = taxonomicGroup,
            quiet = FALSE
        )
        if (is.null(taxonomicGroup)) {
            taxonomicGroup <- basename(dirname(baseUrl))
        }
        if (is.null(genomeBuild)) {
            genomeBuild <- currentRefseqGenomeBuild(
                organism = organism,
                taxonomicGroup = taxonomicGroup
            )
        }
        releaseUrl <- pasteUrl(baseUrl, "all_assembly_versions", genomeBuild)
        outputDir <- initDir(outputDir)
        outputBasename <- kebabCase(tolower(paste(
            organism,
            genomeBuild,
            "refseq",
            release
        )))
        outputDir <- file.path(outputDir, outputBasename)
        h1(sprintf(
            paste(
                "Downloading RefSeq genome for {.emph %s}",
                "%s %d from {.url %s} to {.path %s}."
            ),
            organism,
            genomeBuild,
            release,
            releaseUrl,
            outputDir
        ))
        assert(
            !isADir(outputDir),
            msg = sprintf("Genome exists at '%s'.", outputDir)
        )
        outputDir <- initDir(outputDir)
        args <- list(
            genomeBuild = genomeBuild,
            outputDir = outputDir,
            releaseUrl = releaseUrl,
            cache = cache
        )
        info <- list()
        info[["date"]] <- Sys.Date()
        info[["metadata"]] <-
            do.call(what = .downloadRefSeqMetadata, args = args)
        info[["genome"]] <-
            do.call(what = .downloadRefseqGenome, args = args)
        info[["transcriptome"]] <-
            do.call(what = .downloadRefSeqTranscriptome, args = args)
        info[["annotation"]] <-
            do.call(what = .downloadRefSeqAnnotation, args = args)
        ## Optionally download the no-alt analysis set (ALT-free), useful for
        ## STAR and other aligners that don't handle ALT contigs well.
        if (isTRUE(analysisSet)) {
            ## nolint start
            info[["analysisSet"]] <- .downloadRefseqAnalysisSet(
                organism = organism,
                genomeBuild = genomeBuild,
                outputDir = outputDir,
                cache = cache
            )
            ## nolint end
        }
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
            "RefSeq genome downloaded successfully to {.path %s}.",
            outputDir
        ))
        invisible(info)
    }


## Updated 2022-05-24.
.downloadRefSeqAnnotation <-
    function(genomeBuild, outputDir, releaseUrl, cache) {
        urls <- c(
            gff = pasteUrl(
                releaseUrl,
                paste0(genomeBuild, "_genomic.gff.gz")
            ),
            gtf = pasteUrl(
                releaseUrl,
                paste0(genomeBuild, "_genomic.gtf.gz")
            )
        )
        files <- .downloadUrls(
            urls = urls,
            outputDir = file.path(outputDir, "annotation"),
            cache = cache
        )
        gffFile <- files[["gff"]]
        gtfFile <- files[["gtf"]]
        ## Create relative path symlinks.
        if (!isWindows() && requireNamespace("withr", quietly = TRUE)) {
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
            gffSymlink <- "annotation.gff3.gz"
            gtfSymlink <- "annotation.gtf.gz"
            withr::with_dir(
                new = outputDir,
                code = {
                    file.symlink(from = gffRelativeFile, to = gffSymlink)
                    file.symlink(from = gtfRelativeFile, to = gtfSymlink)
                }
            )
            files[["gffSymlink"]] <- gffSymlink
            files[["gtfSymlink"]] <- gtfSymlink
        }
        ## Save genomic ranges.
        genes <- makeGRangesFromGff(gffFile, level = "genes")
        transcripts <- makeGRangesFromGff(gffFile, level = "transcripts")
        saveRDS(
            object = genes,
            file = file.path(outputDir, "genes.rds")
        )
        saveRDS(
            object = transcripts,
            file = file.path(outputDir, "transcripts.rds")
        )
        ## Save transcript-to-gene mappings.
        t2g <- TxToGene(makeGRangesFromGff(
            file = gffFile,
            level = "transcripts",
            ignoreVersion = FALSE,
            extraMcols = FALSE
        ))
        saveRDS(object = t2g, file = file.path(outputDir, "tx2gene.rds"))
        t2gFile <- export(
            object = t2g,
            con = file.path(outputDir, "tx2gene.csv.gz")
        )
        files[["tx2gene"]] <- t2gFile
        invisible(list(files = files, urls = urls))
    }


## Updated 2022-05-24.
.downloadRefseqGenome <-
    function(genomeBuild, outputDir, releaseUrl, cache) {
        urls <- c(
            fasta = pasteUrl(
                releaseUrl,
                paste0(genomeBuild, "_genomic.fna.gz")
            )
        )
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
            fastaSymlink <- "genome.fa.gz"
            withr::with_dir(
                new = outputDir,
                code = file.symlink(from = fastaRelativeFile, to = fastaSymlink)
            )
            files[["fastaSymlink"]] <- fastaSymlink
        }
        invisible(list(files = files, urls = urls))
    }


## Updated 2023-04-15.
.downloadRefSeqMetadata <-
    function(genomeBuild, outputDir, releaseUrl, cache) {
        urls <- c(
            ## > readme = pasteUrl(releaseUrl, "README.txt"),
            annotationHashes = pasteUrl(releaseUrl, "annotation_hashes.txt"),
            assemblyStatus = pasteUrl(releaseUrl, "assembly_status.txt"),
            md5checksums = pasteUrl(releaseUrl, "md5checksums.txt"),
            assemblyRegions = pasteUrl(
                releaseUrl,
                paste0(genomeBuild, "_assembly_regions.txt")
            ),
            assemblyReport = pasteUrl(
                releaseUrl,
                paste0(genomeBuild, "_assembly_report.txt")
            ),
            assemblyStats = pasteUrl(
                releaseUrl,
                paste0(genomeBuild, "_assembly_stats.txt")
            )
        )
        files <- .downloadUrls(
            urls = urls,
            outputDir = file.path(outputDir, "metadata"),
            cache = cache
        )
        invisible(list(files = files, urls = urls))
    }


## Updated 2022-05-24.
.downloadRefSeqTranscriptome <-
    function(genomeBuild, outputDir, releaseUrl, cache) {
        urls <- c(
            fasta = pasteUrl(
                releaseUrl,
                paste0(genomeBuild, "_rna.fna.gz")
            )
        )
        files <- .downloadUrls(
            urls = urls,
            outputDir = file.path(outputDir, "transcriptome"),
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
            fastaSymlink <- "transcriptome.fa.gz"
            withr::with_dir(
                new = outputDir,
                code = file.symlink(from = fastaRelativeFile, to = fastaSymlink)
            )
            files[["fastaSymlink"]] <- fastaSymlink
        }
        invisible(list(files = files, urls = urls))
    }


## Updated 2026-05-31.
## NCBI alignment pipeline sets: ALT-contig-free FASTA files for Homo sapiens
## and Mus musculus. Useful for STAR alignment without ALT handling.
## See: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/
##      GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/
.downloadRefseqAnalysisSet <-
    function(organism, genomeBuild, outputDir, cache) {
        assert(
            isOrganism(organism),
            isString(genomeBuild),
            isString(outputDir),
            isFlag(cache)
        )
        ## Map known genome builds to their alignment pipeline FASTA URLs.
        ## These are GCA (GenBank) accessions stored under a separate path.
        analysisSetUrls <- list(
            GCF_000001405.40_GRCh38.p14 = pasteUrl(
                "ftp.ncbi.nlm.nih.gov",
                "genomes",
                "all",
                "GCA",
                "000",
                "001",
                "405",
                "GCA_000001405.15_GRCh38",
                "seqs_for_alignment_pipelines.ucsc_ids",
                "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz",
                protocol = "https"
            ),
            GCF_000001635.27_GRCm39 = pasteUrl(
                "ftp.ncbi.nlm.nih.gov",
                "genomes",
                "all",
                "GCA",
                "000",
                "001",
                "635",
                "GCA_000001635.9_GRCm39",
                "seqs_for_alignment_pipelines.ucsc_ids",
                "GCA_000001635.9_GRCm39_no_alt_analysis_set.fna.gz",
                protocol = "https"
            )
        )
        buildPrefix <- sub(
            pattern = "^(GCF_[0-9]+\\.[0-9]+_[^_]+).*$",
            replacement = "\\1",
            x = genomeBuild
        )
        url <- analysisSetUrls[[buildPrefix]]
        if (is.null(url)) {
            alertWarning(sprintf(
                paste0(
                    "No alignment pipeline analysis set available for genome ",
                    "build {.val %s}. Skipping."
                ),
                genomeBuild
            ))
            return(invisible(NULL))
        }
        alert(sprintf(
            "Downloading no-alt analysis set FASTA from {.url %s}.",
            url
        ))
        analysisSetDir <- file.path(outputDir, "genome_analysis_set")
        urls <- c(noAltFasta = url)
        files <- .downloadUrls(
            urls = urls,
            outputDir = analysisSetDir,
            cache = cache
        )
        if (!isWindows() && requireNamespace("withr", quietly = TRUE)) {
            fastaFile <- files[["noAltFasta"]]
            fastaRelativeFile <- sub(
                pattern = paste0("^", outputDir, "/"),
                replacement = "",
                x = fastaFile
            )
            noAltSymlink <- "genome.no_alt.fa.gz"
            withr::with_dir(
                new = outputDir,
                code = file.symlink(from = fastaRelativeFile, to = noAltSymlink)
            )
            files[["fastaSymlink"]] <- noAltSymlink
        }
        invisible(list(files = files, urls = urls))
    }
