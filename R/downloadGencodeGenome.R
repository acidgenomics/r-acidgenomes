#' Download GENCODE reference genome
#'
#' @export
#' @note Updated 2021-08-03.
#'
#' @inheritParams downloadEnsemblGenome
#'
#' @return Invisible `list`.
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
        outputDir = ".",
        cache = TRUE
    ) {
        assert(
            isOrganism(organism),
            isString(genomeBuild, nullOK = TRUE),
            isScalar(release) || is.null(release),
            isString(outputDir),
            isFlag(cache)
        )
        organism <- match.arg(
            arg = organism,
            choices = c("Homo sapiens", "Mus musculus")
        )
        if (is.null(genomeBuild)) {
            genomeBuild <- currentGencodeGenomeBuild(organism)
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
        if (identical(genomeBuild, "GRCh37")) {
            releaseURL <- pasteURL(releaseURL, "GRCh37_mapping")
        }
        outputBasename <- kebabCase(tolower(paste(
            organism, genomeBuild, "gencode", release
        )))
        outputDir <- file.path(outputDir, outputBasename)
        h1(sprintf(
            paste(
                "Downloading GENCODE genome for {.emph %s}",
                " %s %s from {.url %s} to {.path %s}."
            ),
            organism, genomeBuild, as.character(release),
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
            do.call(what = .downloadGencodeMetadata, args = args)
        info[["genome"]] <-
            do.call(what = .downloadGencodeGenome, args = args)
        args <- append(x = args, values = list("release" = release))
        info[["transcriptome"]] <-
            do.call(what = .downloadGencodeTranscriptome, args = args)
        info[["annotation"]] <-
            do.call(what = .downloadGencodeAnnotation, args = args)
        info[["args"]] <- args
        info[["call"]] <- tryCatch(
            expr = standardizeCall(),
            error = function(e) NULL
        )
        info[["sessionInfo"]] <- sessionInfo()
        saveRDS(object = info, file = file.path(outputDir, "metadata.rds"))
        alertSuccess(sprintf(
            "GENCODE genome downloaded successfully to {.path %s}.",
            outputDir
        ))
        invisible(info)
    }



## Updated 2021-08-03.
.downloadGencodeAnnotation <-
    function(
        genomeBuild,
        outputDir,
        release,
        releaseURL,
        cache
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
            outputDir = file.path(outputDir, "annotation"),
            cache = cache
        )
        gffFile <- files[["gff"]]
        gtfFile <- files[["gtf"]]
        ## Create symlinks.
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
            gffSymlink <- paste0("annotation.", fileExt(gffFile))
            gtfSymlink <- paste0("annotation.", fileExt(gtfFile))
            file.symlink(from = gffRelativeFile, to = gffSymlink)
            file.symlink(from = gtfRelativeFile, to = gtfSymlink)
            files[["gffSymlink"]] <- gffSymlink
            files[["gtfSymlink"]] <- gtfSymlink
            setwd(wd)
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
        invisible(list("files" = files, "urls" = urls))
    }



## Updated 2021-08-03.
.downloadGencodeGenome <-
    function(
        genomeBuild,
        outputDir,
        releaseURL,
        cache
    ) {
        urls <- c(
            "fasta" = pasteURL(
                releaseURL,
                paste0(genomeBuild, ".primary_assembly.genome.fa.gz")
            )
        )
        files <- .downloadURLs(
            urls = urls,
            outputDir = file.path(outputDir, "genome"),
            cache = cache
        )
        ## Create symlink.
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
            fastaSymlink <- paste0("genome.", fileExt(fastaFile))
            file.symlink(from = fastaRelativeFile, to = fastaSymlink)
            files[["fastaSymlink"]] <- fastaSymlink
            setwd(wd)
        }
        invisible(list("files" = files, "urls" = urls))
    }



## Updated 2021-08-03.
.downloadGencodeMetadata <-
    function(
        genomeBuild,
        outputDir,
        releaseURL,
        cache
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
            outputDir = file.path(outputDir, "metadata"),
            cache = cache
        )
        invisible(list("files" = files, "urls" = urls))
    }



## Updated 2021-08-03.
.downloadGencodeTranscriptome <-
    function(
        genomeBuild,
        outputDir,
        release,
        releaseURL,
        cache
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
            outputDir = file.path(outputDir, "transcriptome"),
            cache = cache
        )
        fastaFile <- files[["fasta"]]
        ## Save transcript-to-gene mappings.
        tx2gene <- makeTx2GeneFromFASTA(fastaFile)
        saveRDS(object = tx2gene, file = file.path(outputDir, "tx2gene.rds"))
        tx2geneFile <- export(
            object = tx2gene,
            file = file.path(outputDir, "tx2gene.csv.gz")
        )
        files[["tx2gene"]] <- tx2geneFile
        ## Create symlink.
        if (!isWindows()) {
            wd <- getwd()
            setwd(outputDir)
            fastaRelativeFile <- sub(
                pattern = paste0("^", outputDir, "/"),
                replacement = "",
                x = fastaFile
            )
            assert(
                isAFile(fastaFile),
                isAFile(fastaRelativeFile)
            )
            fastaSymlink <- paste0("transcriptome.", fileExt(fastaFile))
            file.symlink(from = fastaRelativeFile, to = fastaSymlink)
            files[["fastaSymlink"]] <- fastaSymlink
            setwd(wd)
        }
        invisible(list("files" = files, "urls" = urls))
    }
