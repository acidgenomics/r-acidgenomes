#' Download Ensembl reference genome
#'
#' @export
#' @note Updated 2021-01-06.
#'
#' @inheritParams AcidRoxygen::params
#' @inheritParams params
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
        genomeBuild,
        release = NULL,
        type = c("all", "transcriptome", "genome", "none"),
        annotation = c("all", "gtf", "gff", "none"),
        outputDir = "."
    ) {
        assert(
            isOrganism(organism),
            isString(genomeBuild),
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
        ## FIXME RETHINK THIS STEP.
        organism <- gsub(pattern = " ", replacement = "_", x = organism)
        baseURL <- "ftp://ftp.ensembl.org/pub"
        if (genomeBuild == "GRCh37") {
            assert(is.null(release))
            baseURL <- pasteURL(baseURL, "grch37", protocol = "none")
            release <- 87L
        }
        if (is.null(release)) {
            release <- currentEnsemblVersion()
        }
        releaseURL <- pasteURL(
            baseURL,
            paste0("release-", release),
            protocol = "none"
        )
        outputBasename = kebabCase(tolower(paste(
            organism, genomeBuild, "ensembl", release
        )))
        outputDir <- file.path(outputDir, outputBasename)
        assert(!isADir(outputDir))
        outputDir <- initDir(outputDir)
        args <- list(
            organism = organism,
            genomeBuild = genomeBuild,
            releaseURL = releaseURL,
            outputDir = outputDir
        )
        if (isTRUE(dlList[["type"]][["genome"]])) {
            do.call(what = .downloadEnsemblGenome, args = args)
        }
        if (isTRUE(dlList[["type"]][["transcriptome"]])) {
            do.call(what = .downloadEnsemblTranscriptome, args = args)
        }
        args <- c(args, release = release)
        if (isTRUE(dlList[["annotation"]][["gtf"]])) {
            do.call(what = .downloadEnsemblGTF, args = args)
        }
        if (isTRUE(dlList[["annotation"]][["gff"]])) {
            do.call(what = .downloadEnsemblGFF, args = args)
        }
        saveRDS(
            object = sessionInfo(),
            file = file.path(outputDir, "sessionInfo.rds")
        )
        cli_alert_success(sprintf(
            "Ensembl genome downloaded successfully to {.path %s}.",
            outputDir
        ))
        invisible(outputDir)
    }



## Updated 2020-12-19.
.downloadEnsemblGenome <-
    function(
        organism,
        genomeBuild,
        releaseURL,
        outputDir
    ) {
        outputDir <- initDir(file.path(outputDir, "genome"))
        baseURL <- pasteURL(
            releaseURL, "fasta", tolower(organism), "dna",
            protocol = "none"
        )
        urls <- c(
            "readme" = pasteURL(baseURL, "README", protocol = "none"),
            "checksums" = pasteURL(baseURL, "CHECKSUMS", protocol = "none")
        )
        if (isSubset(organism, c("Homo_sapiens", "Mus_musculus"))) {
            assembly <- "primary_assembly"
        } else {
            assembly <- "toplevel"
        }
        urls[["fasta"]] <- pasteURL(
            baseURL,
            paste(
                organism, genomeBuild, "dna", assembly, "fa.gz",
                sep = "."
            ),
            protocol = "none"
        )
        destfiles <- vapply(
            X = urls,
            FUN = function(url) {
                file.path(outputDir, basename(url))
            },
            FUN.VALUE = character(1L)
        )
        mapply(
            url = urls,
            destfile = destfiles,
            FUN = download,
            SIMPLIFY = FALSE,
            USE.NAMES = FALSE
        )
        invisible(outputDir)
    }



## Updated 2020-12-19.
.downloadEnsemblTranscriptome <-
    function(
        organism,
        genomeBuild,
        releaseURL,
        outputDir
    ) {
        outputDir <- initDir(file.path(outputDir, "transcriptome"))
        baseURL <- pasteURL(
            releaseURL, "fasta", tolower(organism),
            protocol = "none"
        )
        urls <- character()
        destfiles <- character()
        ## cDNA FASTA.
        cdnaBaseURL <- pasteURL(baseURL, "cdna", protocol = "none")
        cdnaOutputDir <- initDir(file.path(outputDir, "cdna"))
        urls[["cdnaReadme"]] <-
            pasteURL(cdnaBaseURL, "README", protocol = "none")
        destfiles[["cdnaReadme"]] <-
            file.path(cdnaOutputDir, basename(urls[["cdnaReadme"]]))
        urls[["cdnaChecksums"]] <-
            pasteURL(cdnaBaseURL, "CHECKSUMS", protocol = "none")
        destfiles[["cdnaChecksums"]] <-
            file.path(cdnaOutputDir, basename(urls[["cdnaChecksums"]]))
        urls[["cdnaFasta"]] <-
            pasteURL(
                cdnaBaseURL,
                paste(organism, genomeBuild, "cdna.all.fa.gz", sep = "."),
                protocol = "none"
            )
        destfiles[["cdnaFasta"]] <-
            file.path(cdnaOutputDir, basename(urls[["cdnaFasta"]]))
        ## ncRNA FASTA.
        ncrnaBaseURL <- pasteURL(baseURL, "ncrna", protocol = "none")
        ncrnaOutputDir <- initDir(file.path(outputDir, "ncrna"))
        urls[["ncrnaReadme"]] <-
            pasteURL(ncrnaBaseURL, "README", protocol = "none")
        destfiles[["ncrnaReadme"]] <-
            file.path(ncrnaOutputDir, basename(urls[["ncrnaReadme"]]))
        urls[["ncrnaChecksums"]] <-
            pasteURL(ncrnaBaseURL, "CHECKSUMS", protocol = "none")
        destfiles[["ncrnaChecksums"]] <-
            file.path(ncrnaOutputDir, basename(urls[["ncrnaChecksums"]]))
        urls[["ncrnaFasta"]] <-
            pasteURL(
                ncrnaBaseURL,
                paste(organism, genomeBuild, "ncrna.fa.gz", sep = "."),
                protocol = "none"
            )
        destfiles[["ncrnaFasta"]] <-
            file.path(ncrnaOutputDir, basename(urls[["ncrnaFasta"]]))
        assert(identical(names(urls), names(destfiles)))
        mapply(
            url = urls,
            destfile = destfiles,
            FUN = download,
            SIMPLIFY = FALSE,
            USE.NAMES = FALSE
        )
        ## Create a merged transcriptome FASTA.
        cli_alert("Creating a merged transcriptome FASTA file.")
        fastaList <- lapply(
            X = c(
                destfiles[["cdnaFasta"]],
                destfiles[["ncrnaFasta"]]
            ),
            FUN = import,
            format = "lines"
        )
        mergeFasta <- do.call(what = c, args = fastaList)
        mergeFastaFile <- export(
            object = mergeFasta,
            file = file.path(outputDir, "transcriptome.fa.gz"),
            overwrite = TRUE
        )
        makeTx2GeneFileFromFASTA(file = mergeFastaFile, source = "ensembl")
        invisible(outputDir)
    }




## if organism in ("Homo_sapiens", "Mus_musculus"):
##gtf_patch_url = paste_url(
##    base_url,
##    organism
##    + "."
##    + build
##    + "."
##    + release
##    + ".chr_patch_hapl_scaff.gtf.gz",
##)
##download(
##    url=gtf_patch_url, output_dir=output_dir, decompress=decompress
##)

## ftp://ftp.ensembl.org/pub/release-102/gtf/homo_sapiens/Homo_sapiens.GRCh38.102.chr_patch_hapl_scaff.gtf.gz

## ftp://ftp.ensembl.org/pub/release-102/gff3/homo_sapiens/Homo_sapiens.GRCh38.102.chr_patch_hapl_scaff.gff3.gz


## Updated 2021-01-06.
.downloadEnsemblGTF <-
    function(
        organism,
        genomeBuild,
        release,
        releaseURL,
        outputDir
    ) {
        outputDir <- initDir(file.path(outputDir, "gtf"))
        baseURL <- pasteURL(
            releaseURL, "gtf", tolower(organism),
            protocol = "none"
        )
        readmeURL <- pasteURL(baseURL, "README", protocol = "none")
        readmeFile <- file.path(outputDir, basename(readmeURL))
        checksumsURL <- pasteURL(baseURL, "CHECKSUMS", protocol = "none")
        checksumsFile <- file.path(outputDir, basename(checksumsURL))
        gtfURL <- pasteURL(
            baseURL,
            paste(organism, genomeBuild, release, "gtf.gz", sep = "."),
            protocol = "none"
        )
        gtfFile <- file.path(outputDir, basename(gtfURL))
        mapply(
            url = c(
                readmeURL,
                checksumsURL,
                gtfURL
            ),
            destfile = c(
                readmeFile,
                checksumsFile,
                gtfFile
            ),
            FUN = download,
            SIMPLIFY = FALSE,
            USE.NAMES = FALSE
        )
        invisible(outputDir)
    }



## Updated 2021-01-05.
.downloadEnsemblGFF <-
    function(
        organism,
        genomeBuild,
        release,
        releaseURL,
        outputDir
    ) {
        outputDir <- initDir(file.path(outputDir, "gff"))
        baseURL <- pasteURL(
            releaseURL, "gff3", tolower(organism),
            protocol = "none"
        )
        readmeURL <- pasteURL(baseURL, "README", protocol = "none")
        readmeFile <- file.path(outputDir, basename(readmeURL))
        checksumsURL <- pasteURL(baseURL, "CHECKSUMS", protocol = "none")
        checksumsFile <- file.path(outputDir, basename(checksumsURL))
        gffURL <- pasteURL(
            baseURL,
            paste(organism, genomeBuild, release, "gff3.gz", sep = "."),
            protocol = "none"
        )
        gffFile <- file.path(outputDir, basename(gffURL))
        mapply(
            url = c(
                readmeURL,
                checksumsURL,
                gffURL
            ),
            destfile = c(
                readmeFile,
                checksumsFile,
                gffFile
            ),
            FUN = download,
            SIMPLIFY = FALSE,
            USE.NAMES = FALSE
        )
        invisible(outputDir)
    }
