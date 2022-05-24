#' Download GENCODE reference genome
#'
#' @export
#' @note Updated 2022-05-24.
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
    function(organism,
             genomeBuild = NULL,
             release = NULL,
             outputDir = getwd(),
             cache = FALSE) {
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
                "%s %s from {.url %s} to {.path %s}."
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
        info[["genome"]] <-
            do.call(what = .downloadGencodeGenome, args = args)
        args <- append(x = args, values = list("release" = release))
        info[["metadata"]] <-
            do.call(what = .downloadGencodeMetadata, args = args)
        info[["transcriptome"]] <-
            do.call(what = .downloadGencodeTranscriptome, args = args)
        info[["annotation"]] <-
            do.call(
                what = .downloadGencodeAnnotation,
                args = append(
                    x = args,
                    values = list("metadataFiles" = info[["metadata"]])
                )
            )
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



## Updated 2022-05-24.
.downloadGencodeAnnotation <-
    function(genomeBuild,
             metadataFiles,
             outputDir,
             release,
             releaseURL,
             cache) {
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
        assert(
            isAFile(gffFile),
            isAFile(gtfFile)
        )
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
            gffSymlink <- paste0("annotation.", fileExt(gffFile))
            gtfSymlink <- paste0("annotation.", fileExt(gtfFile))
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
        genes <- makeGRangesFromGFF(
            gtfFile,
            level = "genes",
            ignoreVersion = FALSE
        )
        transcripts <- makeGRangesFromGFF(
            gtfFile,
            level = "transcripts",
            ignoreVersion = FALSE
        )
        ## Get Entrez and RefSeq identifier mappings.
        entrezGene <- import(
            file = metadataFiles[["files"]][["entrezGene"]],
            format = "tsv",
            colnames = c("txId", "entrezId")
        )
        entrezGene <- as(entrezGene, "DataFrame")
        entrezGene <- leftJoin(
            x = entrezGene,
            y = mcols(transcripts)[, c("txId", "geneId")],
            by = "txId"
        )
        refSeq <- import(
            file = metadataFiles[["files"]][["refSeq"]],
            format = "tsv",
            colnames = c(
                "txId",
                "refSeqRnaId",
                "refSeqProteinId"
            )
        )
        refSeq <- as(refSeq, "DataFrame")
        refSeq <- leftJoin(
            x = refSeq,
            y = mcols(transcripts)[, c("txId", "geneId")],
            by = "txId"
        )
        ## Add Entrez and RefSeq identifiers to gene metadata.
        mcols <- mcols(genes)
        mcols <- leftJoin(
            x = mcols,
            y = .nest2(
                object = entrezGene,
                by = "geneId",
                exclude = "txId"
            ),
            by = "geneId"
        )
        mcols <- leftJoin(
            x = mcols,
            y = .nest2(
                object = refSeq,
                by = "geneId",
                exclude = "txId"
            ),
            by = "geneId"
        )
        mcols <- mcols[, sort(colnames(mcols))]
        mcols(genes) <- mcols
        ## Add Entrez and RefSeq identifiers to transcript metadata.
        mcols <- mcols(transcripts)
        mcols <- leftJoin(
            x = mcols,
            y = .nest2(
                object = entrezGene,
                by = "txId",
                exclude = "geneId"
            ),
            by = "txId"
        )
        mcols <- leftJoin(
            x = mcols,
            y = .nest2(
                object = refSeq,
                by = "txId",
                exclude = "geneId"
            ),
            by = "txId"
        )
        mcols <- mcols[, sort(colnames(mcols))]
        mcols(transcripts) <- mcols
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



## Updated 2022-05-24.
.downloadGencodeGenome <-
    function(genomeBuild,
             outputDir,
             releaseURL,
             cache) {
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



## Updated 2022-05-03.
.downloadGencodeMetadata <-
    function(genomeBuild,
             outputDir,
             release,
             releaseURL,
             cache) {
        urls <- c(
            "readme" = pasteURL(
                releaseURL,
                switch(
                    EXPR = genomeBuild,
                    "GRCh37" = "_README_GRCh37_mapping.txt",
                    "_README.TXT"
                )
            ),
            "md5sums" = pasteURL(releaseURL, "MD5SUMS"),
            ## TSV (without colnames) mapping transcripts to Entrez gene IDs.
            "entrezGene" = pasteURL(
                releaseURL,
                paste0("gencode.v", release, ".metadata.EntrezGene.gz")
            ),
            ## TSV (without colnames) mapping transcripts to RefSeq IDs.
            "refSeq" = pasteURL(
                releaseURL,
                paste0("gencode.v", release, ".metadata.RefSeq.gz")
            )
        )
        files <- .downloadURLs(
            urls = urls,
            outputDir = file.path(outputDir, "metadata"),
            cache = cache
        )
        invisible(list("files" = files, "urls" = urls))
    }



## Updated 2022-05-24.
.downloadGencodeTranscriptome <-
    function(genomeBuild,
             outputDir,
             release,
             releaseURL,
             cache) {
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
        tx2gene <- makeTx2GeneFromFASTA(
            file = fastaFile,
            ignoreVersion = FALSE
        )
        saveRDS(object = tx2gene, file = file.path(outputDir, "tx2gene.rds"))
        tx2geneFile <- export(
            object = tx2gene,
            con = file.path(outputDir, "tx2gene.csv.gz")
        )
        files[["tx2gene"]] <- tx2geneFile
        ## Create relative path symlink.
        if (!isWindows() && requireNamespace("withr", quietly = TRUE)) {
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
