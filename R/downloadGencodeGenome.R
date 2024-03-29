#' Download GENCODE reference genome
#'
#' @export
#' @note Updated 2023-11-28.
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
#' ## >     release = 44L,
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
            isString(genomeBuild, nullOk = TRUE),
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
        releaseUrl <- pasteUrl(
            "ftp.ebi.ac.uk",
            "pub",
            "databases",
            "gencode",
            paste("Gencode", organismShort, sep = "_"),
            paste("release", release, sep = "_"),
            protocol = "ftp"
        )
        if (identical(genomeBuild, "GRCh37")) {
            releaseUrl <- pasteUrl(releaseUrl, "GRCh37_mapping")
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
            releaseUrl, outputDir
        ))
        assert(
            !isADir(outputDir),
            msg = sprintf("Genome exists at '%s'.", outputDir)
        )
        outputDir <- initDir(outputDir)
        args <- list(
            "genomeBuild" = genomeBuild,
            "outputDir" = outputDir,
            "releaseUrl" = releaseUrl,
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
            error = function(e) {
                NULL
            }
        )
        info[["sessionInfo"]] <- sessionInfo()
        saveRDS(object = info, file = file.path(outputDir, "metadata.rds"))
        alertSuccess(sprintf(
            "GENCODE genome downloaded successfully to {.path %s}.",
            outputDir
        ))
        invisible(info)
    }



## Updated 2023-11-22.
.downloadGencodeAnnotation <-
    function(genomeBuild,
             metadataFiles,
             outputDir,
             release,
             releaseUrl,
             cache) {
        urls <- c(
            "gff" = pasteUrl(
                releaseUrl,
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
            "gtf" = pasteUrl(
                releaseUrl,
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
        files <- .downloadUrls(
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
        genes <- makeGRangesFromGff(
            file = gtfFile,
            level = "genes",
            ignoreVersion = FALSE
        )
        transcripts <- makeGRangesFromGff(
            file = gtfFile,
            level = "transcripts",
            ignoreVersion = FALSE
        )
        t2g <- decode(mcols(transcripts)[, c("txId", "geneId")])
        ## Get NCBI gene and RefSeq identifier mappings.
        ncbiGene <- import(
            con = metadataFiles[["files"]][["ncbiGene"]],
            format = "tsv",
            colnames = c("txId", "ncbiGeneId")
        )
        ncbiGene <- as(ncbiGene, "DFrame")
        ncbiGene <- leftJoin(x = ncbiGene, y = t2g, by = "txId")
        refseq <- import(
            con = metadataFiles[["files"]][["refseq"]],
            format = "tsv",
            colnames = c("txId", "refseqRnaId", "refseqProteinId")
        )
        refseq <- as(refseq, "DFrame")
        refseq <- leftJoin(x = refseq, y = t2g, by = "txId")
        ## Add NCBI and RefSeq identifiers to gene metadata.
        mcols <- decode(mcols(genes))
        if (isSubset("ncbiGeneId", names(mcols))) {
            mcols[["ncbiGeneId"]] <- NULL
        }
        if (isSubset("refseqRnaId", names(mcols))) {
            mcols[["refseqProteinId"]] <- NULL
            mcols[["refseqRnaId"]] <- NULL
        }
        mcols <- leftJoin(
            x = mcols,
            y = .nest2(object = ncbiGene, by = "geneId", exclude = "txId"),
            by = "geneId"
        )
        mcols <- leftJoin(
            x = mcols,
            y = .nest2(object = refseq, by = "geneId", exclude = "txId"),
            by = "geneId"
        )
        mcols <- mcols[, sort(colnames(mcols))]
        mcols(genes) <- encode(mcols)
        ## Add NCBI and RefSeq identifiers to transcript metadata.
        mcols <- decode(mcols(transcripts))
        if (isSubset("ncbiGeneId", names(mcols))) {
            mcols[["ncbiGeneId"]] <- NULL
        }
        if (isSubset("refseqRnaId", names(mcols))) {
            mcols[["refseqProteinId"]] <- NULL
            mcols[["refseqRnaId"]] <- NULL
        }
        mcols <- leftJoin(
            x = mcols,
            y = .nest2(object = ncbiGene, by = "txId", exclude = "geneId"),
            by = "txId"
        )
        mcols <- leftJoin(
            x = mcols,
            y = .nest2(object = refseq, by = "txId", exclude = "geneId"),
            by = "txId"
        )
        mcols <- mcols[, sort(colnames(mcols))]
        mcols(transcripts) <- encode(mcols)
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



## Updated 2023-07-28.
.downloadGencodeGenome <-
    function(genomeBuild,
             outputDir,
             releaseUrl,
             cache) {
        urls <- c(
            "fasta" = pasteUrl(
                releaseUrl,
                paste0(genomeBuild, ".primary_assembly.genome.fa.gz")
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



## Updated 2023-11-28.
.downloadGencodeMetadata <-
    function(genomeBuild,
             outputDir,
             release,
             releaseUrl,
             cache) {
        urls <- c(
            "readme" = pasteUrl(
                releaseUrl,
                switch(
                    EXPR = genomeBuild,
                    "GRCh37" = "_README_GRCh37_mapping.txt",
                    "_README.TXT"
                )
            ),
            "md5sums" = pasteUrl(releaseUrl, "MD5SUMS"),
            ## TSV mapping transcripts to NCBI (Entrez) genes.
            "ncbiGene" = pasteUrl(
                releaseUrl,
                paste0(
                    "gencode.v",
                    release,
                    switch(
                        EXPR = genomeBuild,
                        "GRCh37" = "lift37",
                        ""
                    ),
                    ".metadata.EntrezGene.gz"
                )
            ),
            ## TSV (without colnames) mapping transcripts to RefSeq IDs.
            "refseq" = pasteUrl(
                releaseUrl,
                paste0(
                    "gencode.v",
                    release,
                    switch(
                        EXPR = genomeBuild,
                        "GRCh37" = "lift37",
                        ""
                    ),
                    ".metadata.RefSeq.gz"
                )
            )
        )
        files <- .downloadUrls(
            urls = urls,
            outputDir = file.path(outputDir, "metadata"),
            cache = cache
        )
        invisible(list("files" = files, "urls" = urls))
    }



## Updated 2023-07-28.
##
## Regarding pipe delimiter handling:
## https://github.com/nf-core/rnaseq/issues/864
.downloadGencodeTranscriptome <-
    function(genomeBuild,
             outputDir,
             release,
             releaseUrl,
             cache) {
        urls <- c(
            "fasta" = pasteUrl(
                releaseUrl,
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
        files <- .downloadUrls(
            urls = urls,
            outputDir = file.path(outputDir, "transcriptome"),
            cache = cache
        )
        fastaFile <- files[["fasta"]]
        ## Prepare a fixed FASTA file, without "|" in headers.
        fastaFixedFile <- sub(
            pattern = ".transcripts.",
            replacement = ".transcripts_fixed.",
            x = fastaFile,
            fixed = TRUE
        )
        alert(sprintf(
            "Preparing fixed FASTA file {.file %s} without pipe delimiter.",
            basename(fastaFixedFile)
        ))
        lines <- import(con = fastaFile, format = "lines")
        lines <- sub(
            pattern = "^>([^|]+)\\|.+$",
            replacement = ">\\1",
            x = lines
        )
        export(object = lines, con = fastaFixedFile)
        files[["fastaFixed"]] <- fastaFixedFile
        ## Save transcript-to-gene mappings.
        t2g <- makeTxToGeneFromFasta(
            file = fastaFile,
            ignoreVersion = FALSE
        )
        saveRDS(object = t2g, file = file.path(outputDir, "tx2gene.rds"))
        t2gFile <- export(
            object = t2g,
            con = file.path(outputDir, "tx2gene.csv.gz")
        )
        files[["tx2gene"]] <- t2gFile
        ## Create relative path symlink.
        if (!isWindows() && requireNamespace("withr", quietly = TRUE)) {
            fastaRelativeFile <- sub(
                pattern = paste0("^", outputDir, "/"),
                replacement = "",
                x = fastaFixedFile
            )
            fastaSymlink <- paste0("transcriptome.", fileExt(fastaFixedFile))
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
