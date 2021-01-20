## FIXME CREATE STANDARDIZED SYMLINKS AT TOP LEVEL.
## FIXME DOWNLOAD ANNOTATIONS INTO "ANNOTATIONS" DIR.



#' Download Ensembl reference genome
#'
#' @export
#' @note Updated 2021-01-20.
#'
#' @inheritParams AcidRoxygen::params
#' @inheritParams params
#'
#' @param type `character(1)`.
#'   Genome type to download.
#' @param annotation `character(1)`.
#'   Annotation format to download.
#' @param outputDir `character(1)`.
#'   Output directory path.
#'
#' @examples
#' ## This example is bandwidth intensive.
#' ## > downloadEnsemblGenome(
#' ## >     organism = "Homo sapiens",
#' ## >     genomeBuild = "GRCh38",
#' ## >     release = 100L,
#' ## >     type = "transcriptome",
#' ## >     annotation = "gtf"
#' ## > )
downloadEnsemblGenome <-
    function(
        organism,
        genomeBuild = NULL,
        release = NULL,
        type = c("all", "transcriptome", "genome", "none"),
        annotation = c("all", "gtf", "gff", "none"),
        outputDir = "."
    ) {
        assert(
            isOrganism(organism),
            isString(genomeBuild, nullOK = TRUE),
            isInt(release, nullOK = TRUE),
            isString(outputDir)
        )
        outputDir <- initDir(outputDir)
        type <- match.arg(type)
        annotation <- match.arg(annotation)
        if (type == "none" && annotation == "none") {
            stop("'type' or 'annotation' argument is required.")
        }
        dlList <- list(
            "type" = c(
                "genome" = FALSE,
                "transcriptome" = FALSE
            ),
            "annotation" = c(
                "gff" = FALSE,
                "gtf" = FALSE
            )
        )
        switch(
            EXPR = type,
            "all" = {
                dlList[["type"]][["genome"]] <- TRUE
                dlList[["type"]][["transcriptome"]] <- TRUE
            },
            "genome" = {
                dlList[["type"]][["genome"]] <- TRUE
            },
            "transcriptome" = {
                dlList[["type"]][["transcriptome"]] <- TRUE
            }
        )
        switch(
            EXPR = annotation,
            "all" = {
                dlList[["annotation"]][["gff"]] <- TRUE
                dlList[["annotation"]][["gtf"]] <- TRUE
            },
            "gff" = {
                dlList[["annotation"]][["gff"]] <- TRUE
            },
            "gtf" = {
                dlList[["annotation"]][["gtf"]] <- TRUE
            }
        )
        baseURL <- "ftp://ftp.ensembl.org/pub"
        if (is.null(genomeBuild)) {
            genomeBuild <- currentEnsemblBuild(organism)
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
            "organism" = organism,
            "genomeBuild" = genomeBuild,
            "releaseURL" = releaseURL,
            "outputDir" = outputDir
        )
        out <- list()
        if (isTRUE(dlList[["type"]][["genome"]])) {
            out[["type"]][["genome"]] <-
                do.call(what = .downloadEnsemblGenomeFASTA, args = args)
        }
        if (isTRUE(dlList[["type"]][["transcriptome"]])) {
            out[["type"]][["transcriptome"]] <-
                do.call(what = .downloadEnsemblTranscriptomeFASTA, args = args)
        }
        args <- append(x = args, values = list("release" = release))
        if (isTRUE(dlList[["annotation"]][["gtf"]])) {
            out[["annotation"]][["gtf"]] <-
                do.call(what = .downloadEnsemblGTF, args = args)
        }
        if (isTRUE(dlList[["annotation"]][["gff"]])) {
            out[["annotation"]][["gff"]] <-
                do.call(what = .downloadEnsemblGFF, args = args)
        }
        out[["args"]] <- args
        out[["call"]] <- match.call()
        saveRDS(
            object = sessionInfo(),
            file = file.path(outputDir, "sessionInfo.rds")
        )
        alertSuccess(sprintf(
            "Ensembl genome downloaded successfully to {.path %s}.",
            outputDir
        ))
        invisible(out)
    }



## Updated 2021-01-20.
.downloadEnsemblGFF <-
    function(
        organism,
        genomeBuild,
        release,
        releaseURL,
        outputDir
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
            outputDir = file.path(outputDir, "annotation")
        )
        invisible(files)
    }



## Updated 2021-01-20.
.downloadEnsemblGTF <-
    function(
        organism,
        genomeBuild,
        release,
        releaseURL,
        outputDir
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
            outputDir = file.path(outputDir, "annotation")
        )
        invisible(files)
    }



## Updated 2021-01-20.
.downloadEnsemblGenomeFASTA <-
    function(
        organism,
        genomeBuild,
        releaseURL,
        outputDir
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
        fastaFile <- files[["fasta"]]
        assert(isAFile(fastaFile))
        fastaSymlink <-
            file.path(outputDir, paste0("genome.", fileExt(fastaFile)))
        file.symlink(from = fastaFile, to = fastaSymlink)
        files[["fastaSymlink"]] <- fastaSymlink
        invisible(files)
    }



## Updated 2021-01-20.
.downloadEnsemblTranscriptomeFASTA <-
    function(
        organism,
        genomeBuild,
        releaseURL,
        outputDir
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
        fastaSymlink <- file.path(outputDir, basename(mergeFastaFile))
        file.symlink(from = mergeFastaFile, to = fastaSymlink)
        tx2geneFile <-
            makeTx2GeneFileFromFASTA(file = mergeFastaFile, source = "ensembl")
        tx2geneSymlink <- file.path(outputDir, basename(tx2geneFile))
        file.symlink(from = tx2geneFile, to = tx2geneSymlink)
        files <- list(
            "fasta" = list(
                "cdna" = cdnaFiles,
                "ncrna" = ncrnaFiles,
                "merge" = mergeFastaFile
            ),
            "fastaSymlink" = fastaSymlink,
            "tx2gene" = tx2geneFile,
            "tx2geneSymlink" = tx2geneSymlink
        )
        invisible(files)
    }



#' Download multiple genome files in a single call
#'
#' @note Updated 2021-01-07.
#' @noRd
#'
#' @return `character`
#'   Destination files.
.downloadURLs <- function(urls, outputDir) {
    assert(
        allAreURLs(urls),
        isString(outputDir)
    )
    outputDir <- initDir(outputDir)
    destfiles <- vapply(
        X = urls,
        FUN = function(url) {
            file.path(outputDir, basename(url))
        },
        FUN.VALUE = character(1L)
    )
    assert(identical(names(urls), names(destfiles)))
    out <- mapply(
        url = urls,
        destfile = destfiles,
        FUN = download,
        SIMPLIFY = TRUE,
        USE.NAMES = FALSE
    )
    names(out) <- names(urls)
    assert(allAreFiles(out))
    invisible(out)
}
