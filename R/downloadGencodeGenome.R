#' Download GENCODE reference genome
#'
#' @export
#' @note Updated 2021-01-20.
#'
#' @inheritParams downloadEnsemblGenome
#'
#' @examples
#' ## This example is bandwidth intensive.
#' ## > downloadGencodeGenome(
#' ## >     organism = "Homo sapiens",
#' ## >     genomeBuild = "GRCh38",
#' ## >     release = 36L,
#' ## >     type = "transcriptome",
#' ## >     annotation = "gtf"
#' ## > )
downloadGencodeGenome <-
    function(
        organism,
        genomeBuild = NULL,
        release = NULL,
        outputDir = "."
    ) {
        assert(
            isOrganism(organism),
            isString(genomeBuild, nullOK = TRUE),
            isInt(release, nullOK = TRUE),
            isString(outputDir)
        )
        organism <- match.arg(
            arg = organism,
            choices = c("Homo sapiens", "Mus musculus")
        )
        if (is.null(genomeBuild)) {
            genomeBuild <- currentGencodeBuild(organism)
            genomeBuild <- .simpleGenomeBuild(genomeBuild)
        }
        if (is.null(release)) {
            release <- currentGencodeVersion(organism = organism)
        }
        outputDir <- initDir(outputDir)
        organismShort <- switch(
            EXPR = organism,
            "Homo sapiens" = "human",
            "Mus musculus" = "mouse"
        )
        releaseURL <- pasteURL(
            "ftp://ftp.ebi.ac.uk",
            "pub",
            "databases",
            "gencode",
            paste("Gencode", organismShort, sep = "_"),
            paste("release", release, sep = "_")
        )
        if (genomeBuild == "GRCh37") {
            releaseURL <-
                pasteURL(releaseURL, "GRCh37_mapping")
        }
        outputBasename <- kebabCase(tolower(paste(
            organism, genomeBuild, "gencode", release
        )))
        outputDir <- file.path(outputDir, outputBasename)
        h1(sprintf(
            paste(
                "Downloading GENCODE genome for {.emph %s}",
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
            do.call(what = .downloadGencodeMetadata, args = args)
        info[["genome"]] <-
            do.call(what = .downloadGencodeGenome, args = args)
        args <- append(x = args, values = list("release" = release))
        info[["transcriptome"]] <-
            do.call(what = .downloadGencodeTranscriptome, args = args)
        info[["annotation"]] <-
            do.call(what = .downloadGencodeAnnotation, args = args)
        info[["args"]] <- args
        info[["call"]] <- match.call()
        info[["sessionInfo"]] <- sessionInfo()
        saveRDS(object = info, file = file.path(outputDir, "metadata.rds"))
        alertSuccess(sprintf(
            "GENCODE genome downloaded successfully to {.path %s}.",
            outputDir
        ))
        invisible(info)
    }



## Updated 2021-01-20.
.downloadGencodeAnnotation <-
    function(
        genomeBuild,
        outputDir,
        release,
        releaseURL
    ) {
        urls <- c(
            "gff" = pasteURL(
                releaseURL,
                paste0(
                    "gencode.v",
                    release,
                    switch(
                        EXPR = genomeBuild,
                        "GRCh37" = "lift37",
                        ""
                    ),
                    ".annotation.gff3.gz"
                )
            ),
            "gtf" = pasteURL(
                releaseURL,
                paste0(
                    "gencode.v",
                    release,
                    switch(
                        EXPR = genomeBuild,
                        "GRCh37" = "lift37",
                        ""
                    ),
                    ".annotation.gtf.gz"
                )
            )
        )
        files <- .downloadURLs(
            urls = urls,
            outputDir = file.path(outputDir, "annotation")
        )
        ## Generate GFF symlink.
        gffFile <- files[["gff"]]
        assert(isAFile(gffFile))
        gffSymlink <- file.path(
            outputDir,
            paste0("annotation.", fileExt(gffFile))
        )
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
        invisible(list("files" = files, "urls" = urls))
    }



## Updated 2021-01-20.
.downloadGencodeGenome <-
    function(
        genomeBuild,
        outputDir,
        releaseURL
    ) {
        urls <- c(
            "fasta" = pasteURL(
                releaseURL,
                paste0(genomeBuild, ".primary_assembly.genome.fa.gz")
            )
        )
        files <- .downloadURLs(
            urls = urls,
            outputDir = file.path(outputDir, "genome")
        )
        fastaFile <- files[["fasta"]]
        assert(isAFile(fastaFile))
        fastaSymlink <- file.path(
            outputDir,
            paste0("genome.", fileExt(fastaFile))
        )
        file.symlink(from = fastaFile, to = fastaSymlink)
        files[["fastaSymlink"]] <- fastaSymlink
        invisible(list("files" = files, "urls" = urls))
    }



## Updated 2021-01-20.
.downloadGencodeMetadata <-
    function(
        genomeBuild,
        outputDir,
        releaseURL
    ) {
        urls <- c(
            "readme" = pasteURL(
                releaseURL,
                switch(
                    EXPR = genomeBuild,
                    "GRCh37" = "_README_GRCh37_mapping.txt",
                    "_README.TXT"
                )
            ),
            "md5sums" = pasteURL(releaseURL, "MD5SUMS")
        )
        files <- .downloadURLs(
            urls = urls,
            outputDir = file.path(outputDir, "metadata")
        )
        invisible(list("files" = files, "urls" = urls))
    }



## Updated 2021-01-20.
.downloadGencodeTranscriptome <-
    function(
        genomeBuild,
        outputDir,
        release,
        releaseURL
    ) {
        urls <- c(
            "fasta" = pasteURL(
                releaseURL,
                paste0(
                    "gencode.v", release,
                    switch(
                        EXPR = genomeBuild,
                        "GRCh37" = "lift37",
                        ""
                    ),
                    ".transcripts.fa.gz"
                )
            )
        )
        files <- .downloadURLs(
            urls = urls,
            outputDir = file.path(outputDir, "transcriptome")
        )
        fastaFile <- files[["fasta"]]
        assert(isAFile(fastaFile))
        fastaSymlink <- file.path(
            outputDir,
            paste0("transcriptome.", fileExt(fastaFile))
        )
        file.symlink(from = fastaFile, to = fastaSymlink)
        files[["fastaSymlink"]] <- fastaSymlink
        tx2geneFile <- makeTx2GeneFileFromFASTA(
            file = fastaFile,
            outputFile = file.path(outputDir, "tx2gene.csv.gz"),
            source = "gencode"
        )
        files[["tx2gene"]] <- tx2geneFile
        invisible(list("files" = files, "urls" = urls))
    }
