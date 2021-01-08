#' Download RefSeq reference genome
#'
#' @export
#' @note Updated 2021-01-08.
#'
#' @inheritParams currentGenomeBuild
#' @inheritParams downloadEnsemblGenome
#'
#' @seealso
#' - [RefSeq Genomes FTP server](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/)
#' - [Human Genome Resources at NCBI](https://www.ncbi.nlm.nih.gov/projects/genome/guide/human/)
#' - [Genomes Download (FTP) FAQ](https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/)
#'
#' @examples
#' ## This example is bandwidth intensive.
#' ## > downloadRefSeqGenome(
#' ## >     organism = "Homo sapiens",
#' ## >     taxonomicGroup = "vertebrate_mammalian"
#' ## >     genomeBuild = "GRCh38",
#' ## >     type = "transcriptome",
#' ## >     annotation = "gtf"
#' ## > )
downloadRefSeqGenome <-
    function(
        organism,
        taxonomicGroup = NULL,
        genomeBuild = NULL,
        type = c("all", "transcriptome", "genome", "none"),
        annotation = c("all", "gtf", "gff", "none"),
        outputDir = "."
    ) {
        assert(
            isOrganism(organism),
            isString(taxonomicGroup, nullOK = TRUE),
            isString(genomeBuild, nullOK = TRUE),
            isString(outputDir)
        )
        baseURL <- .getRefSeqGenomeURL(
            organism = organism,
            taxonomicGroup = taxonomicGroup,
            quiet = FALSE
        )
        taxonomicGroup <- basename(dirname(baseURL))
        if (is.null(genomeBuild)) {
            genomeBuild <- currentRefSeqGenomeBuild(
                organism = organism,
                taxonomicGroup = taxonomicGroup
            )
        }
        release <- currentRefSeqVersion()
        releaseURL <- pasteURL(baseURL, "reference", genomeBuild)
        assert(url.exists(releaseURL))
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
        outputBasename <- kebabCase(tolower(paste(
            organism, genomeBuild, "refseq", release
        )))
        outputDir <- file.path(outputDir, outputBasename)
        h1(sprintf(
            paste(
                "Downloading RefSeq genome for {.emph %s}",
                " %s %d from {.url %s} to {.path %s}."
            ),
            organism, genomeBuild, release,
            releaseURL, outputDir
        ))
        assert(!isADir(outputDir))
        outputDir <- initDir(outputDir)
        out <- list()






        ## Download the metadata files first.
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
        out[["metadata"]] <- .downloadURLs(urls = urls, outputDir = outputDir)
        ## Now ready to download individual genome files.
        args <- list(
            genomeBuild = genomeBuild,
            releaseURL = releaseURL,
            outputDir = outputDir
        )
        if (isTRUE(dlList[["type"]][["genome"]])) {
            out[["type"]][["genome"]] <-
                do.call(what = .downloadGencodeGenomeFASTA, args = args)
        }
        args <- c(args, "release" = release)
        if (isTRUE(dlList[["type"]][["transcriptome"]])) {
            out[["type"]][["transcriptome"]] <-
                do.call(what = .downloadGencodeTranscriptomeFASTA, args = args)
        }
        args <- c(args, release = release)
        if (isTRUE(dlList[["annotation"]][["gtf"]])) {
            out[["annotation"]][["gtf"]] <-
                do.call(what = .downloadGencodeGTF, args = args)
        }
        if (isTRUE(dlList[["annotation"]][["gff"]])) {
            out[["annotation"]][["gff"]] <-
                do.call(what = .downloadGencodeGFF, args = args)
        }
        out[["args"]] <- args
        out[["call"]] <- match.call()
        saveRDS(
            object = sessionInfo(),
            file = file.path(outputDir, "sessionInfo.rds")
        )
        alertSuccess(sprintf(
            "GENCODE genome downloaded successfully to {.path %s}.",
            outputDir
        ))
        invisible(out)
    }



## Updated 2021-01-07.
.downloadGencodeGFF <-
    function(
        genomeBuild,
        release,
        releaseURL,
        outputDir
    ) {
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
            )
        )
        .downloadURLs(urls = urls, outputDir = file.path(outputDir, "gff"))
    }



## Updated 2021-01-07.
.downloadGencodeGTF <-
    function(
        genomeBuild,
        release,
        releaseURL,
        outputDir
    ) {
        urls <- c(
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
        .downloadURLs(urls = urls, outputDir = file.path(outputDir, "gtf"))
    }



## Updated 2021-01-07.
.downloadGencodeGenomeFASTA <-
    function(
        genomeBuild,
        releaseURL,
        outputDir
    ) {
        urls <- c(
            "fasta" = pasteURL(
                releaseURL,
                paste0(genomeBuild, ".primary_assembly.genome.fa.gz")
            )
        )
        .downloadURLs(urls = urls, outputDir = file.path(outputDir, "genome"))
    }



## Updated 2021-01-07.
.downloadGencodeTranscriptomeFASTA <-
    function(
        genomeBuild,
        release,
        releaseURL,
        outputDir
    ) {
        out <- list()
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
        out[["fasta"]] <- .downloadURLs(
            urls = urls,
            outputDir = file.path(outputDir, "transcriptome")
        )
        fastaFile <- out[["fasta"]][["fasta"]]
        assert(isAFile(fastaFile))
        out[["tx2gene"]] <- makeTx2GeneFileFromFASTA(
            file = fastaFile,
            source = "gencode"
        )
        invisible(out)
    }
