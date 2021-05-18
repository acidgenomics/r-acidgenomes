#' Download Ensembl reference genome
#'
#' @export
#' @note Updated 2021-02-17.
#'
#' @inheritParams AcidRoxygen::params
#' @inheritParams params
#'
#' @param outputDir `character(1)`.
#'   Output directory path.
#'
#' @return Invisible `list`.
#'
#' @examples
#' ## This example is bandwidth intensive.
#' ## > downloadEnsemblGenome(
#' ## >     organism = "Homo sapiens",
#' ## >     genomeBuild = "GRCh38",
#' ## >     release = 100L
#' ## > )
downloadEnsemblGenome <-
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
        outputDir <- initDir(outputDir)
        baseURL <- "ftp://ftp.ensembl.org/pub"
        if (is.null(genomeBuild)) {
            genomeBuild <- currentEnsemblGenomeBuild(organism)
            genomeBuild <- .simpleGenomeBuild(genomeBuild)
        }
        if (genomeBuild == "GRCh37") {
            assert(is.null(release))
            baseURL <- pasteURL(baseURL, "grch37")
            release <- 87L
        }
        if (is.null(release)) {
            release <- currentEnsemblVersion()
        }
        releaseURL <- pasteURL(baseURL, paste0("release-", release))
        outputBasename <- kebabCase(tolower(paste(
            organism, genomeBuild, "ensembl", release
        )))
        outputDir <- file.path(outputDir, outputBasename)
        h1(sprintf(
            paste(
                "Downloading Ensembl genome for {.emph %s}",
                " %s %d from {.url %s} to {.path %s}."
            ),
            organism, genomeBuild, release,
            releaseURL, outputDir
        ))
        assert(!isADir(outputDir))
        outputDir <- initDir(outputDir)
        args <- list(
            "genomeBuild" = genomeBuild,
            "organism" = organism,
            "outputDir" = outputDir,
            "releaseURL" = releaseURL
        )
        info <- list()
        info[["date"]] <- Sys.Date()
        info[["genome"]] <-
            do.call(what = .downloadEnsemblGenome, args = args)
        info[["transcriptome"]] <-
            do.call(what = .downloadEnsemblTranscriptome, args = args)
        args <- append(x = args, values = list("release" = release))
        info[["annotation"]][["gff"]] <-
            do.call(what = .downloadEnsemblGFF, args = args)
        info[["annotation"]][["gtf"]] <-
            do.call(what = .downloadEnsemblGTF, args = args)
        info[["args"]] <- args
        info[["call"]] <- tryCatch(
            expr = standardizeCall(),
            error = function(e) NULL
        )
        info[["sessionInfo"]] <- sessionInfo()
        saveRDS(object = info, file = file.path(outputDir, "metadata.rds"))
        alertSuccess(sprintf(
            "Ensembl genome downloaded successfully to {.path %s}.",
            outputDir
        ))
        invisible(info)
    }



## Updated 2021-02-17.
.downloadEnsemblGFF <-
    function(
        genomeBuild,
        organism,
        outputDir,
        releaseURL,
        release
    ) {
        baseURL <- pasteURL(releaseURL, "gff3", snakeCase(organism))
        urls <- c(
            "readme" = pasteURL(baseURL, "README"),
            "checksums" = pasteURL(baseURL, "CHECKSUMS"),
            "gff" = pasteURL(
                baseURL,
                paste(
                    gsub(pattern = " ", replacement = "_", x = organism),
                    genomeBuild, release, "gff3.gz",
                    sep = "."
                )
            )
        )
        if (isSubset(organism, c("Homo sapiens", "Mus musculus"))) {
            urls[["gff2"]] <- pasteURL(
                baseURL,
                paste(
                    gsub(pattern = " ", replacement = "_", x = organism),
                    genomeBuild, release, "chr_patch_hapl_scaff", "gff3.gz",
                    sep = "."
                )
            )
        }
        files <- .downloadURLs(
            urls = urls,
            outputDir = file.path(outputDir, "annotation", "gff3")
        )
        ## Create symlink.
        if (!isWindows()) {
            wd <- getwd()
            setwd(outputDir)
            gffFile <- files[["gff"]]
            gffRelativeFile <- sub(
                pattern = paste0("^", outputDir, "/"),
                replacement = "",
                x = gffFile
            )
            assert(
                isAFile(gffFile),
                isAFile(gffRelativeFile)
            )
            gffSymlink <- paste0("annotation.", fileExt(gffFile))
            file.symlink(from = gffRelativeFile, to = gffSymlink)
            files[["gffSymlink"]] <- gffSymlink
            setwd(wd)
        }
        invisible(list("files" = files, "urls" = urls))
    }



## Updated 2021-02-17.
.downloadEnsemblGTF <-
    function(
        genomeBuild,
        organism,
        outputDir,
        release,
        releaseURL
    ) {
        baseURL <- pasteURL(releaseURL, "gtf", snakeCase(organism))
        urls <- c(
            "readme" = pasteURL(baseURL, "README"),
            "checksums" = pasteURL(baseURL, "CHECKSUMS"),
            "gtf" = pasteURL(
                baseURL,
                paste(
                    gsub(pattern = " ", replacement = "_", x = organism),
                    genomeBuild, release, "gtf.gz",
                    sep = "."
                )
            )
        )
        if (isSubset(organism, c("Homo sapiens", "Mus musculus"))) {
            urls[["gtf2"]] <- pasteURL(
                baseURL,
                paste(
                    gsub(pattern = " ", replacement = "_", x = organism),
                    genomeBuild, release, "chr_patch_hapl_scaff", "gtf.gz",
                    sep = "."
                )
            )
        }
        files <- .downloadURLs(
            urls = urls,
            outputDir = file.path(outputDir, "annotation", "gtf")
        )
        gtfFile <- files[["gtf"]]
        ## Create symlink.
        if (!isWindows()) {
            wd <- getwd()
            setwd(outputDir)
            gtfRelativeFile <- sub(
                pattern = paste0("^", outputDir, "/"),
                replacement = "",
                x = gtfFile
            )
            assert(
                isAFile(gtfFile),
                isAFile(gtfRelativeFile)
            )
            gtfSymlink <- paste0("annotation.", fileExt(gtfFile))
            file.symlink(from = gtfRelativeFile, to = gtfSymlink)
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



## Updated 2021-02-17.
.downloadEnsemblGenome <-
    function(
        genomeBuild,
        organism,
        outputDir,
        releaseURL
    ) {
        baseURL <- pasteURL(releaseURL, "fasta", snakeCase(organism), "dna")
        urls <- c(
            "readme" = pasteURL(baseURL, "README"),
            "checksums" = pasteURL(baseURL, "CHECKSUMS")
        )
        if (isSubset(organism, c("Homo sapiens", "Mus musculus"))) {
            assembly <- "primary_assembly"
        } else {
            assembly <- "toplevel"
        }
        urls[["fasta"]] <- pasteURL(
            baseURL,
            paste(
                gsub(pattern = " ", replacement = "_", x = organism),
                genomeBuild,
                "dna",
                assembly,
                "fa.gz",
                sep = "."
            )
        )
        files <- .downloadURLs(
            urls = urls,
            outputDir = file.path(outputDir, "genome")
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



## Updated 2021-02-17.
.downloadEnsemblTranscriptome <-
    function(
        genomeBuild,
        organism,
        outputDir,
        releaseURL
    ) {
        baseURL <- pasteURL(releaseURL, "fasta", snakeCase(organism))
        ## Download cDNA FASTA files.
        cdnaBaseURL <- pasteURL(baseURL, "cdna")
        urls <- c(
            "readme" = pasteURL(cdnaBaseURL, "README"),
            "checksums" = pasteURL(cdnaBaseURL, "CHECKSUMS"),
            "fasta" = pasteURL(
                cdnaBaseURL,
                paste(
                    gsub(pattern = " ", replacement = "_", x = organism),
                    genomeBuild, "cdna", "all", "fa.gz",
                    sep = "."
                )
            )
        )
        cdnaFiles <- .downloadURLs(
            urls = urls,
            outputDir = file.path(outputDir, "transcriptome", "cdna")
        )
        ## Download ncRNA FASTA files.
        ncrnaBaseURL <- pasteURL(baseURL, "ncrna")
        urls <- c(
            "readme" = pasteURL(ncrnaBaseURL, "README"),
            "checksums" = pasteURL(ncrnaBaseURL, "CHECKSUMS"),
            "fasta" = pasteURL(
                ncrnaBaseURL,
                paste(
                    gsub(pattern = " ", replacement = "_", x = organism),
                    genomeBuild, "ncrna", "fa.gz",
                    sep = "."
                )
            )
        )
        ncrnaFiles <- .downloadURLs(
            urls =  urls,
            outputDir = file.path(outputDir, "transcriptome", "ncrna")
        )
        ## Create a merged transcriptome FASTA.
        alert("Creating a merged transcriptome FASTA file.")
        fastaList <- lapply(
            X = c(
                cdnaFiles[["fasta"]],
                ncrnaFiles[["fasta"]]
            ),
            FUN = import,
            format = "lines"
        )
        mergeFasta <- do.call(what = c, args = fastaList)
        mergeFastaFile <- export(
            object = mergeFasta,
            file = file.path(outputDir, "transcriptome", "transcriptome.fa.gz"),
            overwrite = TRUE
        )
        ## Save transcript-to-gene mappings.
        tx2gene <- makeTx2GeneFromFASTA(mergeFastaFile)
        saveRDS(object = tx2gene, file = file.path(outputDir, "tx2gene.rds"))
        tx2geneFile <- export(
            object = tx2gene,
            file = file.path(outputDir, "tx2gene.csv.gz")
        )
        files <- list(
            "fasta" = list(
                "cdna" = cdnaFiles,
                "ncrna" = ncrnaFiles,
                "merge" = mergeFastaFile
            ),
            "tx2gene" = tx2geneFile
        )
        ## Create FASTA symlink.
        if (!isWindows()) {
            wd <- getwd()
            setwd(outputDir)
            fastaFile <- mergeFastaFile
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
