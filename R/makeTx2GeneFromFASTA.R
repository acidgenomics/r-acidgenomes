#' Make a Tx2Gene object from transcriptome FASTA
#'
#' @export
#' @note Updated 2023-07-31.
#'
#' @details
#' RefSeq transcript FASTA (e.g. "GCF_000001405.39_GRCh38.p13_rna.fna.gz")
#' doesn't contain gene identifiers, and is not supported.
#'
#' @inheritParams AcidRoxygen::params
#' @inheritParams params
#'
#' @return `Tx2Gene`.
#'
#' @examples
#' ## Ensembl ====
#' file <- AcidBase::pasteURL(
#'     "ftp.ensembl.org",
#'     "pub",
#'     "release-102",
#'     "fasta",
#'     "homo_sapiens",
#'     "cdna",
#'     "Homo_sapiens.GRCh38.cdna.all.fa.gz",
#'     protocol = "https"
#' )
#' t2g <- makeTx2GeneFromFASTA(file)
#' print(t2g)
#'
#' ## GENCODE ====
#' ## > file <- AcidBase::pasteURL(
#' ## >     "ftp.ebi.ac.uk",
#' ## >     "pub",
#' ## >     "databases",
#' ## >     "gencode",
#' ## >     "Gencode_human",
#' ## >     "release_32",
#' ## >     "gencode.v32.transcripts.fa.gz",
#' ## >     protocol = "https"
#' ## > )
#' ## > t2g <- makeTx2GeneFromFASTA(file)
#' ## > print(t2g)
#'
#' ## FlyBase ====
#' ## > file <- AcidBase::pasteURL(
#' ## >     "ftp.flybase.net",
#' ## >     "releases",
#' ## >     "FB2019_05",
#' ## >     "dmel_r6.30",
#' ## >     "fasta",
#' ## >     "dmel-all-transcript-r6.30.fasta.gz",
#' ## >     protocol = "https"
#' ## > )
#' ## > t2g <- makeTx2GeneFromFASTA(file)
#' ## > print(t2g)
#'
#' ## WormBase ====
#' ## > file <- AcidBase::pasteURL(
#' ## >     "ftp.wormbase.org",
#' ## >     "pub",
#' ## >     "wormbase",
#' ## >     "releases",
#' ## >     "WS272",
#' ## >     "species",
#' ## >     "c_elegans",
#' ## >     "PRJNA13758",
#' ## >     "c_elegans.PRJNA13758.WS272.mRNA_transcripts.fa.gz",
#' ## >     protocol = "https"
#' ## > )
#' ## > t2g <- makeTx2GeneFromFASTA(file)
#' ## > print(t2g)
makeTx2GeneFromFASTA <-
    function(file,
             ignoreVersion = FALSE) {
        assert(
            isString(file),
            isFlag(ignoreVersion)
        )
        if (isAFile(file)) {
            file <- realpath(file)
        }
        alert(sprintf(
            "Making {.cls %s} from FASTA file ({.file %s}).",
            "Tx2Gene", file
        ))
        x <- import(con = .cacheIt(file), format = "lines")
        x <- grep(pattern = "^>", x = x, value = TRUE)
        if (!hasLength(x)) {
            abort(sprintf("Unsupported FASTA: {.file %s}.", basename(file)))
        }
        x <- substr(x, start = 2L, stop = nchar(x))
        ## Detect the provider of the FASTA.
        head <- head(x, n = 10L)
        if (any(grepl(
            pattern = paste0(
                "^(ENS.*T[0-9]{11}\\.[0-9]+)\\s",
                ".+",
                "gene:(ENS.*G[0-9]{11}\\.[0-9]+)"
            ),
            x = head
        ))) {
            ## Note that Ensembl includes "gene:" key.
            ## e.g. "ENST00000632684.1 cdna chromosome.*"
            provider <- "Ensembl"
        } else if (any(grepl(pattern = "FlyBase", x = head))) {
            provider <- "FlyBase"
        } else if (any(grepl(
            pattern = paste0(
                "^(ENS.*T[0-9]{11}\\.[0-9]+)(_PAR_Y)?\\|",
                "(ENS.*G[0-9]{11}\\.[0-9]+)(_PAR_Y)?\\|"
            ),
            x = head
        ))) {
            ## Note that GENCODE uses pipes to separate.
            ## e.g. "ENST00000456328.2|ENSG00000223972.5|.*".
            provider <- "GENCODE"
        } else if (any(grepl(
            pattern = "\\sgene=(WBGene[0-9]{8})$",
            x = head
        ))) {
            provider <- "WormBase"
        } else {
            abort(sprintf(
                "Unsupported FASTA: {.file %s}.",
                basename(file)
            ))
        }
        alertInfo(sprintf("%s transcriptome detected.", provider))
        switch(
            EXPR = provider,
            "Ensembl" = {
                x <- strsplit(x = x, split = " ", fixed = TRUE)
                x <- lapply(
                    X = x,
                    FUN = function(x) {
                        x[c(1L, 4L)]
                    }
                )
                x <- do.call(what = rbind, args = x)
                x[, 2L] <- gsub(
                    pattern = "^gene:",
                    replacement = "",
                    x = x[, 2L]
                )
                assert(
                    allAreMatchingRegex(
                        x = x[, 1L],
                        pattern = "^ENS.*T[0-9]{11}\\.[0-9]+$"
                    ),
                    allAreMatchingRegex(
                        x = x[, 2L],
                        pattern = "^ENS.*G[0-9]{11}\\.[0-9]+$"
                    )
                )
            },
            "FlyBase" = {
                x <- strsplit(x = x, split = " ", fixed = TRUE)
                x <- lapply(
                    X = x,
                    FUN = function(x) {
                        x[c(1L, 9L)]
                    }
                )
                x <- do.call(what = rbind, args = x)
                x[, 2L] <- gsub(
                    pattern = "^.*\\b(FBgn[0-9]{7})\\b.*$",
                    replacement = "\\1",
                    x = x[, 2L]
                )
                assert(
                    allAreMatchingRegex(
                        x = x[, 1L],
                        pattern = "^FBtr[0-9]{7}$"
                    ),
                    allAreMatchingRegex(
                        x = x[, 2L],
                        pattern = "^FBgn[0-9]{7}$"
                    )
                )
            },
            "GENCODE" = {
                x <- strsplit(x = x, split = "|", fixed = TRUE)
                x <- lapply(
                    X = x,
                    FUN = function(x) {
                        x[c(1L, 2L)]
                    }
                )
                x <- do.call(what = rbind, args = x)
                assert(
                    allAreMatchingRegex(
                        x = x[, 1L],
                        pattern = "^ENS.*T[0-9]{11}\\.[0-9]+(_PAR_Y)?$"
                    ),
                    allAreMatchingRegex(
                        x = x[, 2L],
                        pattern = "^ENS.*G[0-9]{11}\\.[0-9]+(_PAR_Y)?$"
                    )
                )
            },
            "WormBase" = {
                x <- strsplit(x = x, split = " ", fixed = TRUE)
                x <- lapply(
                    X = x,
                    FUN = function(x) {
                        x[c(1L, 2L)]
                    }
                )
                x <- do.call(what = rbind, args = x)
                x[, 2L] <- gsub(
                    pattern = "^gene=",
                    replacement = "",
                    x = x[, 2L]
                )
                assert(
                    allAreMatchingRegex(
                        x = x[, 1L],
                        pattern = "^[A-Za-z0-9_]+\\.[a-z0-9]+\\.[0-9]+$"
                    ),
                    allAreMatchingRegex(
                        x = x[, 2L],
                        pattern = "^WBGene[0-9]{8}$"
                    )
                )
            },
            abort(sprintf("Unsupported FASTA: {.file %s}.", basename(file)))
        )
        if (isTRUE(ignoreVersion)) {
            x[, 1L] <- stripTranscriptVersions(x[, 1L])
            x[, 2L] <- stripGeneVersions(x[, 2L])
        }
        out <- Tx2Gene(x)
        meta <- list(
            "call" = tryCatch(
                expr = standardizeCall(),
                error = function(e) {
                    NULL
                }
            ),
            "date" = Sys.Date(),
            "file" = file,
            "ignoreVersion" = ignoreVersion,
            "packageVersion" = .pkgVersion,
            "provider" = provider
        )
        metadata(out) <- meta
        out
    }
