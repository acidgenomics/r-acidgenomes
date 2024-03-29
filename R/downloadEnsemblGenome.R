#' Download Ensembl reference genome
#'
#' @export
#' @note Updated 2023-10-12.
#'
#' @inheritParams AcidRoxygen::params
#' @inheritParams params
#'
#' @param outputDir `character(1)`.
#' Output directory path.
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
    function(organism,
             genomeBuild = NULL,
             release = NULL,
             outputDir = getwd(),
             cache = FALSE) {
        assert(
            isOrganism(organism),
            isString(genomeBuild, nullOk = TRUE),
            isInt(release, nullOk = TRUE),
            isString(outputDir),
            isFlag(cache)
        )
        outputDir <- initDir(outputDir)
        ## Can switch to HTTPS if FTP server is down, but often slower.
        baseUrl <- pasteUrl("ftp.ensembl.org", "pub", protocol = "ftp")
        if (is.null(genomeBuild)) {
            genomeBuild <- currentEnsemblGenomeBuild(organism)
            genomeBuild <- .simpleGenomeBuild(genomeBuild)
        }
        if (identical(genomeBuild, "GRCh37")) {
            assert(is.null(release))
            baseUrl <- pasteUrl(baseUrl, "grch37")
            release <- 87L
        }
        if (is.null(release)) {
            release <- currentEnsemblVersion()
        }
        releaseUrl <- pasteUrl(baseUrl, paste0("release-", release))
        outputBasename <- kebabCase(tolower(paste(
            organism, genomeBuild, "ensembl", release
        )))
        outputDir <- file.path(outputDir, outputBasename)
        h1(sprintf(
            paste(
                "Downloading Ensembl genome for {.emph %s}",
                "%s %d from {.url %s} to {.path %s}."
            ),
            organism, genomeBuild, release,
            releaseUrl, outputDir
        ))
        assert(
            !isADir(outputDir),
            msg = sprintf("Genome exists at {.path %s}.", outputDir)
        )
        outputDir <- initDir(outputDir)
        args <- list(
            "organism" = organism,
            "genomeBuild" = genomeBuild,
            "outputDir" = outputDir,
            "releaseUrl" = releaseUrl,
            "cache" = cache
        )
        info <- list()
        info[["date"]] <- Sys.Date()
        info[["genome"]] <-
            do.call(what = .downloadEnsemblGenome, args = args)
        info[["transcriptome"]] <-
            do.call(what = .downloadEnsemblTranscriptome, args = args)
        args <- append(x = args, values = list("release" = release))
        info[["annotation"]][["gff"]] <-
            do.call(what = .downloadEnsemblGff, args = args)
        info[["annotation"]][["gtf"]] <-
            do.call(what = .downloadEnsemblGtf, args = args)
        if (!isSubset(genomeBuild, "GRCh37")) {
            info[["metadata"]] <-
                do.call(what = .downloadEnsemblMetadata, args = args)
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
            "Ensembl genome downloaded successfully to {.path %s}.",
            outputDir
        ))
        invisible(info)
    }



## Updated 2023-04-14.
.downloadEnsemblGff <-
    function(organism,
             genomeBuild,
             release,
             outputDir,
             releaseUrl,
             cache) {
        baseUrl <- pasteUrl(releaseUrl, "gff3", snakeCase(organism))
        urls <- c(
            "readme" = pasteUrl(baseUrl, "README"),
            "checksums" = pasteUrl(baseUrl, "CHECKSUMS"),
            "gff" = pasteUrl(
                baseUrl,
                paste(
                    gsub(pattern = " ", replacement = "_", x = organism),
                    genomeBuild, release, "gff3.gz",
                    sep = "."
                )
            )
        )
        if (
            isSubset(
                x = organism,
                y = c("Homo sapiens", "Mus musculus")
            ) &&
                ## NOTE Not supported for new GRCm39 build.
                isSubset(
                    x = genomeBuild,
                    y = c("GRCh37", "GRCh38", "GRCm38")
                )
        ) {
            urls[["gff2"]] <- pasteUrl(
                baseUrl,
                paste(
                    gsub(pattern = " ", replacement = "_", x = organism),
                    genomeBuild, release, "chr_patch_hapl_scaff", "gff3.gz",
                    sep = "."
                )
            )
        }
        files <- .downloadUrls(
            urls = urls,
            outputDir = file.path(outputDir, "annotation", "gff3"),
            cache = cache
        )
        ## Create relative path symlink.
        if (!isWindows() && requireNamespace("withr", quietly = TRUE)) {
            gffFile <- files[["gff"]]
            gffRelativeFile <- sub(
                pattern = paste0("^", outputDir, "/"),
                replacement = "",
                x = gffFile
            )
            gffSymlink <- paste0("annotation.", fileExt(gffFile))
            withr::with_dir(
                new = outputDir,
                code = {
                    file.symlink(from = gffRelativeFile, to = gffSymlink)
                }
            )
            files[["gffSymlink"]] <- gffSymlink
        }
        invisible(list("files" = files, "urls" = urls))
    }



## Updated 2023-04-14.
.downloadEnsemblGtf <-
    function(organism,
             genomeBuild,
             release,
             outputDir,
             releaseUrl,
             cache) {
        baseUrl <- pasteUrl(releaseUrl, "gtf", snakeCase(organism))
        urls <- c(
            "readme" = pasteUrl(baseUrl, "README"),
            "checksums" = pasteUrl(baseUrl, "CHECKSUMS"),
            "gtf" = pasteUrl(
                baseUrl,
                paste(
                    gsub(pattern = " ", replacement = "_", x = organism),
                    genomeBuild, release, "gtf.gz",
                    sep = "."
                )
            )
        )
        if (
            isSubset(
                x = organism,
                y = c("Homo sapiens", "Mus musculus")
            ) &&
                ## NOTE Not supported for new GRCm39 build.
                isSubset(
                    x = genomeBuild,
                    y = c("GRCh37", "GRCh38", "GRCm38")
                )
        ) {
            urls[["gtf2"]] <- pasteUrl(
                baseUrl,
                paste(
                    gsub(pattern = " ", replacement = "_", x = organism),
                    genomeBuild, release, "chr_patch_hapl_scaff", "gtf.gz",
                    sep = "."
                )
            )
        }
        files <- .downloadUrls(
            urls = urls,
            outputDir = file.path(outputDir, "annotation", "gtf"),
            cache = cache
        )
        gtfFile <- files[["gtf"]]
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



## Updated 2023-04-14.
.downloadEnsemblGenome <-
    function(organism,
             genomeBuild,
             outputDir,
             releaseUrl,
             cache) {
        baseUrl <- pasteUrl(releaseUrl, "fasta", snakeCase(organism), "dna")
        urls <- c(
            "readme" = pasteUrl(baseUrl, "README"),
            "checksums" = pasteUrl(baseUrl, "CHECKSUMS")
        )
        if (isSubset(organism, c("Homo sapiens", "Mus musculus"))) {
            assembly <- "primary_assembly"
        } else {
            assembly <- "toplevel"
        }
        urls[["fasta"]] <- pasteUrl(
            baseUrl,
            paste(
                gsub(pattern = " ", replacement = "_", x = organism),
                genomeBuild,
                "dna",
                assembly,
                "fa.gz",
                sep = "."
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



## Updated 2023-04-14.
.downloadEnsemblMetadata <-
    function(organism,
             genomeBuild,
             release,
             outputDir,
             releaseUrl,
             cache) {
        baseUrl <- pasteUrl(releaseUrl, "tsv", snakeCase(organism))
        organism2 <- gsub(pattern = " ", replacement = "_", x = organism)
        urls <- c(
            "ena" = pasteUrl(
                baseUrl,
                paste(
                    organism2, genomeBuild, release, "ena", "tsv", "gz",
                    sep = "."
                )
            ),
            "entrez" = pasteUrl(
                baseUrl,
                paste(
                    organism2, genomeBuild, release, "entrez", "tsv", "gz",
                    sep = "."
                )
            ),
            "karyotype" = pasteUrl(
                baseUrl,
                paste(
                    organism2, genomeBuild, release, "karyotype", "tsv", "gz",
                    sep = "."
                )
            ),
            "refseq" = pasteUrl(
                baseUrl,
                paste(
                    organism2, genomeBuild, release, "refseq", "tsv", "gz",
                    sep = "."
                )
            ),
            "uniprot" = pasteUrl(
                baseUrl,
                paste(
                    organism2, genomeBuild, release, "uniprot", "tsv", "gz",
                    sep = "."
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



## Updated 2022-05-24.
.downloadEnsemblTranscriptome <-
    function(organism,
             genomeBuild,
             outputDir,
             releaseUrl,
             cache) {
        baseUrl <- pasteUrl(releaseUrl, "fasta", snakeCase(organism))
        ## Download cDNA FASTA files.
        cdnaBaseUrl <- pasteUrl(baseUrl, "cdna")
        urls <- c(
            "readme" = pasteUrl(cdnaBaseUrl, "README"),
            "checksums" = pasteUrl(cdnaBaseUrl, "CHECKSUMS"),
            "fasta" = pasteUrl(
                cdnaBaseUrl,
                paste(
                    gsub(pattern = " ", replacement = "_", x = organism),
                    genomeBuild, "cdna", "all", "fa.gz",
                    sep = "."
                )
            )
        )
        cdnaFiles <- .downloadUrls(
            urls = urls,
            outputDir = file.path(outputDir, "transcriptome", "cdna"),
            cache = cache
        )
        ## Download ncRNA FASTA files.
        ncrnaBaseUrl <- pasteUrl(baseUrl, "ncrna")
        urls <- c(
            "readme" = pasteUrl(ncrnaBaseUrl, "README"),
            "checksums" = pasteUrl(ncrnaBaseUrl, "CHECKSUMS"),
            "fasta" = pasteUrl(
                ncrnaBaseUrl,
                paste(
                    gsub(pattern = " ", replacement = "_", x = organism),
                    genomeBuild, "ncrna", "fa.gz",
                    sep = "."
                )
            )
        )
        ncrnaFiles <- .downloadUrls(
            urls = urls,
            outputDir = file.path(outputDir, "transcriptome", "ncrna"),
            cache = cache
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
            con = file.path(outputDir, "transcriptome", "transcriptome.fa.gz"),
            overwrite = TRUE
        )
        ## Save transcript-to-gene mappings.
        t2g <- makeTxToGeneFromFasta(
            file = mergeFastaFile,
            ignoreVersion = FALSE
        )
        saveRDS(object = t2g, file = file.path(outputDir, "tx2gene.rds"))
        t2gFile <- export(
            object = t2g,
            con = file.path(outputDir, "tx2gene.csv.gz")
        )
        files <- list(
            "fasta" = list(
                "cdna" = cdnaFiles,
                "ncrna" = ncrnaFiles,
                "merge" = mergeFastaFile
            ),
            "tx2gene" = t2gFile
        )
        ## Create relative path symlink.
        if (!isWindows() && requireNamespace("withr", quietly = TRUE)) {
            fastaFile <- mergeFastaFile
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
