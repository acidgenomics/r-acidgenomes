## FIXME THIS ISNT PRODUCTION READY YET.



#' Download GENCODE reference genome
#'
#' @export
#' @note Updated 2021-01-05.
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
        genomeBuild,
        release = NULL,
        type = c("all", "transcriptome", "genome", "none"),
        annotation = c("all", "gtf", "gff", "none"),
        outputDir = "."
    ) {
        assert(
            isString(organism),
            isString(genomeBuild),
            isInt(release, nullOK = TRUE),
            isString(outputDir)
        )
        organism <- match.arg(
            arg = organism,
            choices = c("Homo sapiens", "Mus musculus")
        )
        ## FIXME Consider reworking with `currentGencodeBuild` function.
        genomeBuild <- match.arg(
            arg = genomeBuild,
            choices = switch(
                EXPR = organism,
                "Homo sapiens" = c("GRCh38", "GRCh37"),
                "Mus musculus" = "GRCm38"
            )
        )
        urls <- character()
        destfiles <- character()
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
        organism <- gsub(pattern = " ", replacement = "_", x = organism)
        baseURL <- pasteURL(
            "ftp://ftp.ebi.ac.uk",
            "pub",
            "databases",
            "gencode",
            paste("Gencode", organismShort, sep = "_"),
            paste("release", release, sep = "_"),
            protocol = "none"
        )
        outputBasename <- kebabCase(tolower(paste(
            organism, genomeBuild, "gencode", release
        )))
        outputDir <- file.path(outputDir, outputBasename)
        assert(!isADir(outputDir))
        outputDir <- initDir(outputDir)
        if (genomeBuild == "GRCh37") {
            baseURL <- pasteURL(baseURL, "GRCh37_mapping", protocol = "none")
        }
        ## README file.
        urls[["readme"]] <- pasteURL(
            baseURL,
            switch(
                EXPR = genomeBuild,
                "GRCh37" = "_README_GRCh37_mapping.txt",
                "_README.TXT"
            ),
            protocol = "none"
        )
        destfiles[["readme"]] <- file.path(outputDir, basename(urls[["readme"]]))
        ## MD5 checksums file.
        urls[["md5sums"]] <- pasteURL(baseURL, "MD5SUMS", protocol = "none")
        destfiles[["md5sums"]] <- file.path(
            outputDir, basename(urls[["md5sums"]])
        )
        ## Genome FASTA file.
        if (isTRUE(dlList[["type"]][["genome"]])) {
            urls[["genome"]] <- pasteURL(
                baseURL,
                paste0(genomeBuild, ".primary_assembly.genome.fa.gz"),
                protocol = "none"
            )
            destfiles[["genome"]] <- file.path(
                outputDir, "genome", basename(urls[["genome"]])
            )
        }
        ## Transcriptome FASTA file.
        if (isTRUE(dlList[["type"]][["transcriptome"]])) {
            urls[["transcriptome"]] <- pasteURL(
                baseURL,
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
            destfiles[["transcriptome"]] <- file.path(
                outputDir, "transcriptome", basename(urls[["transcriptome"]])
            )
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

        ## FIXME RETHINK TX2GENE HANDLING HERE.

        saveRDS(
            object = sessionInfo(),
            file = file.path(outputDir, "sessionInfo.rds")
        )
        cli_alert_success(sprintf(
            "GENCODE genome downloaded successfully to {.path %s}.",
            outputDir
        ))
        invisible(outputDir)
    }



## Updated 2021-01-05.
.downloadGencodeGenome <-
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
        readmeURL <- pasteURL(baseURL, "README", protocol = "none")
        readmeFile <- file.path(outputDir, basename(readmeURL))
        checksumsURL <- pasteURL(baseURL, "CHECKSUMS", protocol = "none")
        checksumsFile <- file.path(outputDir, basename(checksumsURL))
        if (isSubset(organism, c("Homo_sapiens", "Mus_musculus"))) {
            assembly <- "primary_assembly"
        } else {
            assembly <- "toplevel"
        }
        fastaURL <- pasteURL(
            baseURL,
            paste(
                organism, genomeBuild, "dna", assembly, "fa.gz",
                sep = "."
            ),
            protocol = "none"
        )
        fastaFile <- file.path(outputDir, basename(fastaURL))
        mapply(
            url = c(
                readmeURL,
                checksumsURL,
                fastaURL
            ),
            destfile = c(
                readmeFile,
                checksumsFile,
                fastaFile
            ),
            FUN = download,
            SIMPLIFY = FALSE,
            USE.NAMES = FALSE
        )
        invisible(outputDir)
    }



## Updated 2021-01-05.
.downloadGencodeTranscriptome <-
    function(
        organism,
        genomeBuild,
        releaseURL,
        outputDir
    ) {
        outputDir = initDir(file.path(outputDir, "transcriptome"))
        baseURL <- pasteURL(
            releaseURL, "fasta", tolower(organism),
            protocol = "none"
        )
        ## cDNA FASTA.
        cdnaOutputDir <- initDir(file.path(outputDir, "cdna"))
        cdnaBaseURL <-
            pasteURL(baseURL, "cdna", protocol = "none")
        cdnaReadmeURL <-
            pasteURL(cdnaBaseURL, "README", protocol = "none")
        cdnaReadmeFile <-
            file.path(cdnaOutputDir, basename(cdnaReadmeURL))
        cdnaChecksumsURL <-
            pasteURL(cdnaBaseURL, "CHECKSUMS", protocol = "none")
        cdnaChecksumsFile <-
            file.path(cdnaOutputDir, basename(cdnaChecksumsURL))
        cdnaFastaURL <-
            pasteURL(
                cdnaBaseURL,
                paste(organism, genomeBuild, "cdna.all.fa.gz", sep = "."),
                protocol = "none"
            )
        cdnaFastaFile <-
            file.path(cdnaOutputDir, basename(cdnaFastaURL))
        ## ncRNA FASTA.
        ncrnaOutputDir <- initDir(file.path(outputDir, "ncrna"))
        ncrnaBaseURL <-
            pasteURL(baseURL, "ncrna", protocol = "none")
        ncrnaReadmeURL <-
            pasteURL(ncrnaBaseURL, "README", protocol = "none")
        ncrnaReadmeFile <-
            file.path(ncrnaOutputDir, basename(ncrnaReadmeURL))
        ncrnaChecksumsURL <-
            pasteURL(ncrnaBaseURL, "CHECKSUMS", protocol = "none")
        ncrnaChecksumsFile <-
            file.path(ncrnaOutputDir, basename(ncrnaChecksumsURL))
        ncrnaFastaURL <- pasteURL(
            ncrnaBaseURL,
            paste(organism, genomeBuild, "ncrna.fa.gz", sep = "."),
            protocol = "none"
        )
        ncrnaFastaFile <-
            file.path(ncrnaOutputDir, basename(ncrnaFastaURL))
        mapply(
            url = c(
                cdnaReadmeURL,
                cdnaChecksumsURL,
                cdnaFastaURL,
                ncrnaReadmeURL,
                ncrnaChecksumsURL,
                ncrnaFastaURL
            ),
            destfile = c(
                cdnaReadmeFile,
                cdnaChecksumsFile,
                cdnaFastaFile,
                ncrnaReadmeFile,
                ncrnaChecksumsFile,
                ncrnaFastaFile
            ),
            FUN = download,
            SIMPLIFY = FALSE,
            USE.NAMES = FALSE
        )
        # Create a merged transcriptome FASTA.
        cdnaFasta <- import(file = cdnaFastaFile, format = "lines")
        ncrnaFasta <- import(file = ncrnaFastaFile, format = "lines")
        mergeFasta <- c(cdnaFasta, ncrnaFasta)
        mergeFastaFile <- export(
            object = mergeFasta,
            file = file.path(outputDir, "transcriptome.fa"),
            overwrite = TRUE
        )
        mergeFastaFile <- compress(
            file = mergeFastaFile,
            ext = "gz",
            remove = !decompress,
            overwrite = TRUE
        )
        tx2gene <- makeTx2GeneFromFASTA(
            file = mergeFastaFile,
            source = "gencode"
        )
        tx2geneFile <- file.path(outputDir, "tx2gene.csv")
        if (isFALSE(decompress)) {
            tx2geneFile <- paste0(tx2geneFile, ".gz")
        }
        export(
            object = tx2gene,
            file = tx2geneFile,
            overwrite = TRUE
        )
        invisible(outputDir)
    }



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
