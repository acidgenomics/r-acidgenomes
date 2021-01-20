## FIXME NEED TO PARSE THE REFSEQ GTF HEADER TO GET THE RELEASE VERSION...
## FIXME GET THE RELEASE VERSION FROM THE GTF METADATA?
##       OTHERWISE THIS CAN BE INCORRECT...
## FIXME NEED TO GET TRANSCRIPTS TO FLAT GRANGES, THEN WE CAN ADD TX2GENE
##       SUPPORT? Alternatively, don't output a tx2gene for this automatically?
## FIXME CREATE STANDARDIZED SYMLINKS AT TOP LEVEL.
## FIXME DOWNLOAD ANNOTATIONS INTO "ANNOTATIONS" DIR.



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
#' @note Updated 2021-01-20.
#'
#' @inheritParams currentGenomeBuild
#' @inheritParams downloadEnsemblGenome
#'
#' @seealso
#' - [Human Genome Resources at NCBI](https://www.ncbi.nlm.nih.gov/projects/genome/guide/human/)
#' - [RefSeq Genomes FTP server](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/)
#' - [Genomes Download (FTP) FAQ](https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/)
#' - [bcbio hg38 reference genome](https://steinbaugh.com/posts/bcbio-hg38.html)
#' - [Heng Li: Which human reference genome to use?](http://lh3.github.io/2017/11/13/which-human-reference-genome-to-use)
#' - [GRCh38 assembly for alignment pipelines](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/)
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
        info[["metadata"]] <-
            do.call(what = .downloadRefSeqMetadata, args = args)
        info[["genome"]] <-
            do.call(what = .downloadRefSeqGenome, args = args)
        info[["transcriptome"]] <-
            do.call(what = .downloadRefSeqTranscriptome, args = args)
        args <- c(args, release = release)
        out[["annotation"]] <-
            do.call(what = .downloadRefSeqAnnotation, args = args)
        ## Export transcript-to-gene mappings.


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
        ## FIXME INCLUDE SYMLINKS.
        ##
        ##
        ## FIXME MAKE THE TX2GENE FROM THE GFF?
        ## > tx2gene <- makeTx2GeneFileFromGFF(file = out[["annotation"]][["gff"]])
        ##
        invisible(list("files" = files, "urls" = urls))
    }

## FIXME MOVE THIS INTO ANNOTATION.
## Updated 2021-01-08.
.downloadRefSeqGTF <-
    function(
        releaseURL,
        genomeBuild,
        outputDir
    ) {
        urls <- c(

        )
        .downloadURLs(urls = urls, outputDir = file.path(outputDir, "gtf"))
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
        ## FIXME Need to add symlink.
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
        ## FIXME NEED TO HANDLE TX2GENE HERE...

        ## FIXME Need to add symlink.
        invisible(list("files" = files, "urls" = urls))
    }
