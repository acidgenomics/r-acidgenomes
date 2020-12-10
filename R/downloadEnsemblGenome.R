#' Download Ensembl reference genome
#'
#' @export
#' @note Updated 2020-12-10.
#'
#' @inheritParams AcidRoxygen::params
#' @inheritParams params
#'
#' @param releaseURL `integer(1)` or `NULL`.
#' @param outputDir `character(1)`.
#' @param decompress `integer(1)`.
#'
#' @examples
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
        ## FIXME Tighten these up by sharing args to function calls.
        if (isTRUE(dlList[["type"]][["genome"]])) {
            .downloadEnsemblGenome(
                organism = organism,
                genomeBuild = genomeBuild,
                releaseURL = releaseURL,
                outputDir = outputDir,
                decompress = decompress
            )
        }
        if (isTRUE(dlList[["type"]][["transcriptome"]])) {
            .downloadEnsemblTranscriptome(
                organism = organism,
                genomeBuild = genomeBuild,
                releaseURL = releaseURL,
                outputDir = outputDir,
                decompress = decompress
            )
        }
        if (isTRUE(dlList[["annotation"]][["gtf"]])) {
            .downloadEnsemblGTF(
                organism = organism,
                genomeBuild = genomeBuild,
                release = release,
                releaseURL = releaseURL,
                outputDir = outputDir,
                decompress = decompress
            )
        }
        if (isTRUE(dlList[["annotation"]][["gff"]])) {
            .downloadEnsemblGFF(
                organism = organism,
                genomeBuild = genomeBuild,
                release = release,
                releaseURL = releaseURL,
                outputDir = outputDir,
                decompress = decompress
            )
        }
        cli_alert_success(sprintf(
            "Ensembl genome downloaded successfully to {.path %s}.",
            outputDir
        ))
        invisible(outputDir)
    }



## Updated 2020-12-10.
.downloadEnsemblGenome <-
    function(
        organism,
        genomeBuild,
        releaseURL,
        outputDir,
        decompress
    ) {
        outputDir = initDir(file.path(outputDir, "genome"))
        baseURL = pasteURL(
            releaseURL, "fasta", tolower(organism), "dna",
            protocol = "none"
        )
        readmeURL = pasteURL(baseURL, "README", protocol = "none")
        checksumsURL = pasteURL(baseURL, "CHECKSUMS", protocol = "none")
        if (isSubset(organism, c("Homo_sapiens", "Mus_musculus"))) {
            assembly <- "primary_assembly"
        } else {
            assembly <- "toplevel"
        }
        fastaURL <- pasteURL(
            baseURL,
            paste0(organism, ".", genomeBuild, ".dna.", assembly, ".fa.gz"),
            protocol = "none"
        )
        readmeFile <- file.path(outputDir, basename(readmeURL))
        checksumsFile <- file.path(outputDir, basename(checksumsURL))
        fastaFile <- file.path(outputDir, basename(fastaURL))
        ## FIXME IMPORT `download` from pipette here.
        ## FIXME DEFINE THIS IN ACIDBASE, AS A HARDENED FUNCTION.
        timeout <- getOption("timeout")
        ## FIXME ONLY SET IF TIMEOUT IS A NUMBER.
        options("timeout" = 99999L)
        status <- mapply(
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
            FUN = download.file,
            SIMPLIFY = FALSE,
            USE.NAMES = FALSE
        )
        assert(all(unlist(status) == 0L))
        if (isTRUE(decompress)) {
            ## FIXME DECOMPRESS THE FASTA, IF DESIRED.
        }
        options("timeout" = timeout)
    }



## Updated 2020-12-10.
downloadEnsemblTranscriptome <-
    function(
        organism,
        genomeBuild,
        releaseURL,
        outputDir,
        decompress
    ) {
        # output_dir = join(output_dir, "transcriptome")
        # transcriptome_file = join(output_dir, "transcriptome.fa.gz")
        # base_url = paste_url(release_url, "fasta", organism.lower())
        # # cDNA FASTA.
        # cdna_output_dir = join(output_dir, "cdna")
        # cdna_base_url = paste_url(base_url, "cdna")
        # cdna_readme_url = paste_url(cdna_base_url, "README")
        # cdna_checksums_url = paste_url(cdna_base_url, "CHECKSUMS")
        # cdna_fasta_file_basename = organism + "." + genomeBuild + ".cdna.all.fa.gz"
        # cdna_fasta_file = join(cdna_output_dir, cdna_fasta_file_basename)
        # cdna_fasta_url = paste_url(cdna_base_url, cdna_fasta_file_basename)
        # download(url=cdna_readme_url, output_dir=cdna_output_dir)
        # download(url=cdna_checksums_url, output_dir=cdna_output_dir)
        # download(
        #     url=cdna_fasta_url, output_file=cdna_fasta_file, decompress=decompress
        # )
        # # ncRNA FASTA.
        # ncrna_output_dir = join(output_dir, "ncrna")
        # ncrna_base_url = paste_url(base_url, "ncrna")
        # ncrna_readme_url = paste_url(ncrna_base_url, "README")
        # ncrna_checksums_url = paste_url(ncrna_base_url, "CHECKSUMS")
        # ncrna_fasta_file_basename = organism + "." + genomeBuild + ".ncrna.fa.gz"
        # ncrna_fasta_file = join(ncrna_output_dir, ncrna_fasta_file_basename)
        # ncrna_fasta_url = paste_url(ncrna_base_url, ncrna_fasta_file_basename)
        # download(url=ncrna_readme_url, output_dir=ncrna_output_dir)
        # download(url=ncrna_checksums_url, output_dir=ncrna_output_dir)
        # download(
        #     url=ncrna_fasta_url,
        #     output_file=ncrna_fasta_file,
        #     decompress=decompress,
        # )
        # # Merged transcriptome FASTA.
        # #
        # # This method is memory efficient. It automatically reads the input files
        # # chunk by chunk for you, which is more more efficient and reading the
        # # input files in and will work even if some of the input files are too
        # # large to fit into memory.
        # #
        # # See also:
        # # - https://stackoverflow.com/a/18209002
        # # - https://stackoverflow.com/a/27077437
        # FIXME REWORK USING CONCATENATE CALL IN R.
        # with open(transcriptome_file, "wb") as output_file:
        #     for file in [cdna_fasta_file, ncrna_fasta_file]:
        #     with open(file, "rb") as file_open:
        #     copyfileobj(file_open, output_file)

        # FIXME GENERATE THE TX2GENE FILE
        makeTx2geneFromFASTA(
            sourceName = "ensembl",
            outputDir = outputDir
        )
    }



## Updated 2020-12-10.
downloadEnsemblGTF <-
    function(
        organism,
        genomeBuild,
        release,
        releaseURL,
        outputDir,
        decompress
    ) {
        # output_dir = join(output_dir, "gtf")
        # base_url = paste_url(release_url, "gtf", organism.lower())
        # readme_url = paste_url(base_url, "README")
        # checksums_url = paste_url(base_url, "CHECKSUMS")
        # gtf_url = paste_url(
        #     base_url, organism + "." + genomeBuild + "." + release + ".gtf.gz"
        # )
        # download(url=readme_url, output_dir=output_dir)
        # download(url=checksums_url, output_dir=output_dir)
        # download(url=gtf_url, output_dir=output_dir, decompress=decompress)
        # if organism in ("Homo_sapiens", "Mus_musculus"):
        #     gtf_patch_url = paste_url(
        #         base_url,
        #         organism
        #         + "."
        #         + genomeBuild
        #         + "."
        #         + release
        #         + ".chr_patch_hapl_scaff.gtf.gz",
        #     )
        # download(
        #     url=gtf_patch_url, output_dir=output_dir, decompress=decompress
        # )
    }



## Updated 2020-12-10.
downloadEnsemblGFF <-
    function(
        organism,
        genomeBuild,
        release,
        releaseURL,
        outputDir,
        decompress
    ) {
        # output_dir = join(output_dir, "gff")
        # base_url = paste_url(release_url, "gff3", organism.lower())
        # readme_url = paste_url(base_url, "README")
        # checksums_url = paste_url(base_url, "CHECKSUMS")
        # gff_url = paste_url(
        #     base_url, organism + "." + genomeBuild + "." + release + ".gff3.gz"
        # )
        # download(url=readme_url, output_dir=output_dir)
        # download(url=checksums_url, output_dir=output_dir)
        # download(url=gff_url, output_dir=output_dir, decompress=decompress)
        # if organism in ("Homo sapiens", "Mus musculus"):
        #     gtf_patch_url = paste_url(
        #         base_url,
        #         organism
        #         + "."
        #         + genomeBuild
        #         + "."
        #         + release
        #         + ".chr_patch_hapl_scaff.gff3.gz",
        #     )
        # download(
        #     url=gtf_patch_url, output_dir=output_dir, decompress=decompress
        # )
    }
