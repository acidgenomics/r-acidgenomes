#' Make a TxToGene object from transcriptome FASTA
#'
#' @export
#' @note Updated 2024-01-02.
#'
#' @details
#' RefSeq transcript FASTA (e.g. "GCF_000001405.39_GRCh38.p13_rna.fna.gz")
#' doesn't contain gene identifiers, and is not supported.
#'
#' @inheritParams AcidRoxygen::params
#' @inheritParams params
#'
#' @return `TxToGene`.
#'
#' @examples
#' ## Ensembl ====
#' file <- AcidBase::pasteUrl(
#'     "ftp.ensembl.org",
#'     "pub",
#'     "release-102",
#'     "fasta",
#'     "homo_sapiens",
#'     "cdna",
#'     "Homo_sapiens.GRCh38.cdna.all.fa.gz",
#'     protocol = "ftp"
#' )
#' t2g <- makeTxToGeneFromFasta(file)
#' print(t2g)
#'
#' ## GENCODE ====
#' ## GRCh38:
#' ## > file <- AcidBase::pasteUrl(
#' ## >     "ftp.ebi.ac.uk",
#' ## >     "pub",
#' ## >     "databases",
#' ## >     "gencode",
#' ## >     "Gencode_human",
#' ## >     "release_32",
#' ## >     "gencode.v32.transcripts.fa.gz",
#' ## >     protocol = "ftp"
#' ## > )
#' ## > t2g <- makeTxToGeneFromFasta(file)
#' ## > print(t2g)
#' ##
#' ## GRCh37:
#' ## > file <- AcidBase::pasteUrl(
#' ## >     "ftp.ebi.ac.uk",
#' ## >     "pub",
#' ## >     "databases",
#' ## >     "gencode",
#' ## >     "Gencode_human",
#' ## >     "release_44",
#' ## >     "GRCh37_mapping",
#' ## >     "gencode.v44lift37.transcripts.fa.gz",
#' ## >     protocol = "ftp"
#' ## > )
#' ## > t2g <- makeTxToGeneFromFasta(file)
#' ## > print(t2g)
#'
#' ## FlyBase ====
#' ## > file <- AcidBase::pasteUrl(
#' ## >     "ftp.flybase.net",
#' ## >     "releases",
#' ## >     "FB2019_05",
#' ## >     "dmel_r6.30",
#' ## >     "fasta",
#' ## >     "dmel-all-transcript-r6.30.fasta.gz",
#' ## >     protocol = "ftp"
#' ## > )
#' ## > t2g <- makeTxToGeneFromFasta(file)
#' ## > print(t2g)
#'
#' ## WormBase ====
#' ## > file <- AcidBase::pasteUrl(
#' ## >     "ftp.wormbase.org",
#' ## >     "pub",
#' ## >     "wormbase",
#' ## >     "releases",
#' ## >     "WS272",
#' ## >     "species",
#' ## >     "c_elegans",
#' ## >     "PRJNA13758",
#' ## >     "c_elegans.PRJNA13758.WS272.mRNA_transcripts.fa.gz",
#' ## >     protocol = "ftp"
#' ## > )
#' ## > t2g <- makeTxToGeneFromFasta(file)
#' ## > print(t2g)
makeTxToGeneFromFasta <-
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
            "TxToGene", file
        ))
        lines <- import(con = .cacheIt(file), format = "lines")
        lines <- grep(pattern = "^>", x = lines, value = TRUE)
        assert(
            hasLength(lines),
            msg = sprintf("Unsupported FASTA: {.file %s}.", basename(file))
        )
        lines <- substr(x = lines, start = 2L, stop = nchar(lines))
        ## FIXME This likely doesn't work for all files but does work for Ensembl.
        first <- strSplit(
            x = lines[[1L]],
            split = " ",
            fixed = TRUE,
            n = 8L
        )[1L, ]
        if (grepl(
            pattern = paste0(
                "^(ENS.*T[0-9]{11}\\.[0-9_]+)(_PAR_Y)?\\|",
                "(ENS.*G[0-9]{11}\\.[0-9_]+)(_PAR_Y)?\\|"
            ),
            x = first[[1L]]
        )) {
            ## GENCODE uses pipes to separate.
            ## e.g. "ENST00000456328.2|ENSG00000223972.5|.*".
            ## Note that GRCh37 identifiers also contain "_".
            ## e.g. "ENST00000456328.2_1".
            provider <- "GENCODE"
        } else if (
            any(grepl(pattern = "FlyBase", x = head)) ||
            any(grepl(pattern = "release=r[.0-9]+; species=Dmel", x = head))
        ) {
            ## FIXME Rework this check with FlyBase files.
            provider <- "FlyBase"
        } else if (any(grepl(
            pattern = "\\sgene=(WBGene[0-9]{8})$",
            x = head
        ))) {
            provider <- "WormBase"
        } else if (
            grepl(
                pattern = " chromosome:",
                x = lines[[1L]],
                fixed = TRUE
            ) &&
            grepl(
                pattern = " gene:",
                x = lines[[1L]],
                fixed = TRUE
            ) &&
            grepl(
                pattern = " gene_biotype:",
                x = lines[[1L]],
                fixed = TRUE
            ) &&
            grepl(
                pattern = " transcript_biotype:",
                x = lines[[1L]],
                fixed = TRUE
            ) &&
            grepl(
                pattern = " gene_symbol:",
                x = lines[[1L]],
                fixed = TRUE
            ) &&
            grepl(
                pattern = " description:",
                x = lines[[1L]],
                fixed = TRUE
            )
        ) {
            ## Alternatively can check for these keys:
            ## - "chromosome:"
            ## - "gene:"
            ## - "gene_biotype:"
            ## - "transcript_biotype:"
            ## - "gene_symbol:"
            ## - "description:"
            ## [1] "ENST00000415118.1"
            ## [2] "cdna"
            ## [3] "chromosome:GRCh38:14:22438547:22438554:1"
            ## [4] "gene:ENSG00000223997.1"
            ## [5] "gene_biotype:TR_D_gene"
            ## [6] "transcript_biotype:TR_D_gene"
            ## [7] "gene_symbol:TRDD1"
            ## [8] "description:T cell receptor delta diversity 1 [Source:HGNC Symbol;Acc:HGNC:12254]"
            provider <- "Ensembl"
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
                ## FIXME Alternatively, we can import as a text connection
                ## delimited by space. This is likely faster, and then we can
                ## just pick columns 1 and 4.
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
                ## FIXME Relax this check, so we can support other genomes,
                ## such as fly and worm.
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
                        pattern = "^ENS.*T[0-9]{11}\\.[0-9_]+(_PAR_Y)?$"
                    ),
                    allAreMatchingRegex(
                        x = x[, 2L],
                        pattern = "^ENS.*G[0-9]{11}\\.[0-9_]+(_PAR_Y)?$"
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
        out <- TxToGene(x)
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
