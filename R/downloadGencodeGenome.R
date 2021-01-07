## FIXME RETHINK THIS...DONT MODIFY
## > organism <- gsub(pattern = " ", replacement = "_", x = organism)



#' Download GENCODE reference genome
#'
#' @export
#' @note Updated 2021-01-07.
#'
#' @inheritParams AcidRoxygen::params
#' @inheritParams params
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
            paste("release", release, sep = "_"),
            protocol = "none"
        )
        if (genomeBuild == "GRCh37") {
            releaseURL <-
                pasteURL(releaseURL, "GRCh37_mapping", protocol = "none")
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
        urls <- c(
            "readme" = pasteURL(
                releaseURL,
                switch(
                    EXPR = genomeBuild,
                    "GRCh37" = "_README_GRCh37_mapping.txt",
                    "_README.TXT"
                ),
                protocol = "none"
            ),
            "md5sums" = pasteURL(releaseURL, "MD5SUMS", protocol = "none")
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



        ## FIXME ===============================================================
        ## THESE NEED TO MIGRATE TO SPECIFIC FUNCTIONS BELOW.

        ## Transcriptome FASTA file.
        ## FIXME NEED TO REWORK TX2GENE PROCESSING HERE.
        if (isTRUE(dlList[["type"]][["transcriptome"]])) {

        }
        ## GTF file.
        if (isTRUE(dlList[["annotation"]][["gtf"]])) {
            urls[["gtf"]] <- pasteURL(
                baseURL,
                paste0(
                    "gencode.v",
                    release,
                    switch(
                        EXPR = genomeBuild,
                        "GRCh37" = "lift37",
                        ""
                    ),
                    ".annotation.gtf.gz"
                ),
                protocol = "none"
            )
            destfiles[["gtf"]] <- file.path(
                outputDir, "gtf", basename(urls[["gtf"]])
            )
        }
        ## GFF3 file.
        if (isTRUE(dlList[["annotation"]][["gff"]])) {
            urls[["gff"]] <- pasteURL(
                baseURL,
                paste0(
                    "gencode.v",
                    release,
                    switch(
                        EXPR = genomeBuild,
                        "GRCh37" = "lift37",
                        ""
                    ),
                    ".annotation.gff3.gz"
                ),
                protocol = "none"
            )
            destfiles[["gff"]] <- file.path(
                outputDir, "gff", basename(urls[["gff"]])
            )
        }
        stopifnot(identical(names(urls), names(destfiles)))
        mapply(
            url = urls,
            destfile = destfiles,
            FUN = download,
            SIMPLIFY = FALSE,
            USE.NAMES = FALSE
        )
        ## FIXME ===============================================================



        args <- list(
            organism = organism,
            genomeBuild = genomeBuild,
            releaseURL = releaseURL,
            outputDir = outputDir
        )
        if (isTRUE(dlList[["type"]][["genome"]])) {
            do.call(what = .downloadGencodeGenome, args = args)
        }

        ## FIXME THESE NEED RELEASE VERSION?

        if (isTRUE(dlList[["type"]][["transcriptome"]])) {
            do.call(what = .downloadGencodeTranscriptome, args = args)
        }
        args <- c(args, release = release)
        if (isTRUE(dlList[["annotation"]][["gtf"]])) {
            do.call(what = .downloadGencodeGTF, args = args)
        }
        if (isTRUE(dlList[["annotation"]][["gff"]])) {
            do.call(what = .downloadGencodeGFF, args = args)
        }
        saveRDS(
            object = sessionInfo(),
            file = file.path(outputDir, "sessionInfo.rds")
        )
        alertSuccess(sprintf(
            "GENCODE genome downloaded successfully to {.path %s}.",
            outputDir
        ))
        invisible(outputDir)
    }



## Updated 2021-01-07.
.downloadGencodeGenome <-
    function(
        organism,
        genomeBuild,
        releaseURL,
        outputDir
    ) {
        outputDir <- initDir(file.path(outputDir, "genome"))
        urls <- c(
            "genome" = pasteURL(
                releaseURL,
                paste0(genomeBuild, ".primary_assembly.genome.fa.gz"),
                protocol = "none"
            )
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



## Updated 2021-01-07.
.downloadGencodeTranscriptome <-
    function(
        organism,
        genomeBuild,
        releaseURL,
        outputDir
    ) {
        outputDir = initDir(file.path(outputDir, "transcriptome"))
        urls <- c(
            "transcriptome" = pasteURL(
                releaseURL,
                paste0(
                    "gencode.v", release,
                    switch(
                        EXPR = genomeBuild,
                        "GRCh37" = "lift37",
                        ""
                    ),
                    ".transcripts.fa.gz"
                ),
                protocol = "none"
            )
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



## FIXME
## Updated 2021-01-05.
.downloadGencodeGTF <-
    function(
        organism,
        genomeBuild,
        release,
        releaseURL,
        outputDir,
        decompress
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
        if (isTRUE(decompress)) {
            decompress(file = gtfFile, remove = FALSE, overwrite = TRUE)
        }
        invisible(outputDir)
    }



## FIXME
## Updated 2021-01-05.
.downloadGencodeGFF <-
    function(
        organism,
        genomeBuild,
        release,
        releaseURL,
        outputDir,
        decompress
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
        if (isTRUE(decompress)) {
            decompress(file = gffFile, remove = FALSE, overwrite = TRUE)
        }
        invisible(outputDir)
    }
