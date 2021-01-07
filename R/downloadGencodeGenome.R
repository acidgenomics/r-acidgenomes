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
        ## Download the metadata files first.
        ## FIXME MAKE THIS A SHARED FUNCTION.
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
        destfiles <- vapply(
            X = urls,
            FUN = function(url) {
                file.path(outputDir, basename(url))
            },
            FUN.VALUE = character(1L)
        )
        ## FIXME MAKE THIS A SHARED FUNCTION.
        mapply(
            url = urls,
            destfile = destfiles,
            FUN = download,
            SIMPLIFY = FALSE,
            USE.NAMES = FALSE
        )
        ## Now ready to download individual genome files.
        args <- list(
            genomeBuild = genomeBuild,
            releaseURL = releaseURL,
            outputDir = outputDir
        )
        if (isTRUE(dlList[["type"]][["genome"]])) {
            do.call(what = .downloadGencodeGenome, args = args)
        }
        args <- c(args, "release" = release)
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
        genomeBuild,
        releaseURL,
        outputDir
    ) {
        outputDir <- initDir(file.path(outputDir, "genome"))
        urls <- c(
            "genome" = pasteURL(
                releaseURL,
                paste0(genomeBuild, ".primary_assembly.genome.fa.gz")
            )
        )
        destfiles <- vapply(
            X = urls,
            FUN = function(url) {
                file.path(outputDir, basename(url))
            },
            FUN.VALUE = character(1L)
        )
        ## FIXME MAKE THIS A SHARED FUNCTION.
        mapply(
            url = urls,
            destfile = destfiles,
            FUN = download,
            SIMPLIFY = FALSE,
            USE.NAMES = FALSE
        )
        invisible(outputDir)
    }



## FIXME NEED TO GENERATE TX2GENE HERE AUTOMATICALLY.
## Updated 2021-01-07.
.downloadGencodeTranscriptome <-
    function(
        genomeBuild,
        release,
        releaseURL,
        outputDir
    ) {
        outputDir <- initDir(file.path(outputDir, "transcriptome"))
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
                )
            )
        )
        destfiles <- vapply(
            X = urls,
            FUN = function(url) {
                file.path(outputDir, basename(url))
            },
            FUN.VALUE = character(1L)
        )
        ## FIXME MAKE THIS A SHARED FUNCTION.
        mapply(
            url = urls,
            destfile = destfiles,
            FUN = download,
            SIMPLIFY = FALSE,
            USE.NAMES = FALSE
        )
        makeTx2GeneFileFromFASTA(
            file = destfiles[["transcriptome"]],
            source = "gencode"
        )
        invisible(outputDir)
    }



## Updated 2021-01-07.
.downloadGencodeGTF <-
    function(
        genomeBuild,
        release,
        releaseURL,
        outputDir
    ) {
        outputDir <- initDir(file.path(outputDir, "gtf"))
        urls <- c(
            "gtf" = pasteURL(
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
                )
            )
        )
        destfiles <- vapply(
            X = urls,
            FUN = function(url) {
                file.path(outputDir, basename(url))
            },
            FUN.VALUE = character(1L)
        )
        ## FIXME MAKE THIS A SHARED FUNCTION.
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
## Updated 2021-01-07.
.downloadGencodeGFF <-
    function(
        genomeBuild,
        release,
        releaseURL,
        outputDir
    ) {
        outputDir <- initDir(file.path(outputDir, "gff"))
        ## FIXME MAKE THIS A SHARED FUNCTION.
        urls <- c(
            "gff" = pasteURL(
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
                )
            )
        )
        destfiles <- vapply(
            X = urls,
            FUN = function(url) {
                file.path(outputDir, basename(url))
            },
            FUN.VALUE = character(1L)
        )
        ## FIXME MAKE THIS A SHARED FUNCTION.
        mapply(
            url = urls,
            destfile = destfiles,
            FUN = download,
            SIMPLIFY = FALSE,
            USE.NAMES = FALSE
        )
        invisible(outputDir)
    }
