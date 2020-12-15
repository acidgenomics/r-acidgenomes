#' Download GENCODE reference genome
#'
#' @export
#' @note Updated 2020-12-14.
#'
#' @inheritParams AcidRoxygen::params
#' @inheritParams params
#'
#' @param releaseURL `integer(1)` or `NULL`.
#' @param outputDir `character(1)`.
#' @param decompress `integer(1)`.
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
downloadGencodeGenome <-
    function(
        organism,
        genomeBuild,
        release = NULL,
        type = c("transcriptome", "genome", "all", "none"),
        annotation = c("gtf", "gff", "all", "none"),
        outputDir = ".",
        decompress = FALSE
    ) {
        assert(
            isString(organism),
            isString(genomeBuild),
            isInt(release, nullOK = TRUE),
            isString(outputDir),
            isFlag(decompress)
        )
        organism <- match.arg(
            arg = organism,
            choices = c("Homo sapiens", "Mus musculus")
        )
        genomeBuild <- match.arg(
            arg = genomeBuild,
            choices = switch(
                EXPR = organism,
                "Homo sapiens" = c("GRCh38", "GRCh37"),
                "Mus musculus" = "GRCm38"
            )
        )
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
            transcriptomeFastaURL <- pasteURL(
                baseURL,
                paste0("gencode.v", release, "lift37.transcripts.fa.gz"),
                protocol = "none"
            )
            gtfURL <- pasteURL(
                baseURL,
                paste0("gencode.v", release, "lift37.annotation.gtf.gz"),
                protocol = "none"
            )
            gffURL <- pasteURL(
                baseURL,
                paste0("gencode.v", release, "lift37.annotation.gff3.gz"),
                protocol = "none"
            )
            readmeURL <- pasteURL(
                baseURL,
                "_README_GRCh37_mapping.txt",
                protocol = "none"
            )
        } else {
            transcriptomeFastaURL <- pasteURL(
                baseURL,
                paste0("gencode.v", release, ".transcripts.fa.gz"),
                protocol = "none"
            )
            gtfURL <- pasteURL(
                baseURL,
                paste0("gencode.v", release, ".annotation.gtf.gz"),
                protocol = "none"
            )
            gffURL <- pasteURL(
                baseURL,
                paste0("gencode.v", release, ".annotation.gff3.gz"),
                protocol = "none"
            )
            readmeURL <- pasteURL(baseURL, "_README.TXT", protocol = "none")
        }
        genomeFastaURL <- pasteURL(
            baseURL,
            paste0(build, ".primary_assembly.genome.fa.gz"),
            protocol = "none"
        )
        md5sumsURL <- pasteURL(baseURL, "MD5SUMS", protocol = "none")
        download(
            url = readmeURL,
            destfile = file.path(outputDir, basename(md5sumsURL))
        )
        download(
            url = md5sumsURL,
            destfile = file.path(outputDir, basename(md5sumsURL))
        )




        if (isTRUE(dlList[["type"]][["genome"]])) {
            ## FIXME In this function, make sure target directory exists
            ## automatically.
            download(
                url = genome_fasta_url,
                output_dir = join(output_dir, "genome"),
                decompress = decompress,
            )
        }
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
        cli_alert_success(sprintf(
            "GENCODE genome downloaded successfully to {.path %s}.",
            outputDir
        ))
        invisible(outputDir)
    }



## Updated 2020-12-15.
.downloadGencodeGenome <-
    function(
        organism,
        genomeBuild,
        releaseURL,
        outputDir,
        decompress
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
        if (isTRUE(decompress)) {
            decompress(file = fastaFile, remove = FALSE, overwrite = TRUE)
        }
        invisible(outputDir)
    }



## Updated 2020-12-10.
.downloadEnsemblTranscriptome <-
    function(
        organism,
        genomeBuild,
        releaseURL,
        outputDir,
        decompress
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
            source = "ensembl"
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



## Updated 2020-12-10.
.downloadEnsemblGTF <-
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



## Updated 2020-12-10.
.downloadEnsemblGFF <-
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
