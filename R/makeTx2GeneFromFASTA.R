#' Make a Tx2Gene object from transcriptome FASTA
#'
#' @export
#' @note RefSeq transcript FASTA
#'   (e.g. "GCF_000001405.39_GRCh38.p13_rna.fna.gz") doesn't contain gene
#'   identifiers, and is not supported.
#' @note Updated 2021-01-29.
#'
#' @inheritParams AcidRoxygen::params
#'
#' @return `Tx2Gene`.
#'
#' @examples
#' ## Ensembl ====
#' file <- pasteURL(
#'     "ftp.ensembl.org",
#'     "pub",
#'     "release-102",
#'     "fasta",
#'     "homo_sapiens",
#'     "cdna",
#'     "Homo_sapiens.GRCh38.cdna.all.fa.gz",
#'     protocol = "ftp"
#' )
#' t2g <- makeTx2GeneFromFASTA(file)
#' print(t2g)
#'
#' ## GENCODE ====
#' ## > file <- pasteURL(
#' ## >     "ftp.ebi.ac.uk",
#' ## >     "pub",
#' ## >     "databases",
#' ## >     "gencode",
#' ## >     "Gencode_human",
#' ## >     "release_32",
#' ## >     "gencode.v32.transcripts.fa.gz",
#' ## >     protocol = "ftp"
#' ## > )
#' ## > t2g <- makeTx2GeneFromFASTA(file)
#' ## > print(t2g)
#'
#' ## FlyBase ====
#' ## > file <- pasteURL(
#' ## >     "ftp.flybase.net",
#' ## >     "releases",
#' ## >     "FB2019_05",
#' ## >     "dmel_r6.30",
#' ## >     "fasta",
#' ## >     "dmel-all-transcript-r6.30.fasta.gz",
#' ## >     protocol = "ftp"
#' ## > )
#' ## > t2g <- makeTx2GeneFromFASTA(file)
#' ## > print(t2g)
#'
#' ## WormBase ====
#' ## > file <- pasteURL(
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
#' ## > t2g <- makeTx2GeneFromFASTA(file)
#' ## > print(t2g)
makeTx2GeneFromFASTA <- function(file) {
    alert(sprintf(
        "Making {.var %s} from FASTA file ({.file %s}).",
        "Tx2Gene", basename(file)
    ))
    x <- import(file = .cacheIt(file), format = "lines")
    x <- grep(pattern = "^>", x = x, value = TRUE)
    if (!hasLength(x)) {
        stop(sprintf("Unsupported FASTA: '%s'.", basename(file)))
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
        stop(sprintf("Unsupported FASTA: '%s'.", basename(file)))
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
            x[, 2L] <- gsub(pattern = "^gene:", replacement = "", x = x[, 2L])
            assert(
                allAreMatchingRegex(
                    pattern = "^ENS.*T[0-9]{11}\\.[0-9]+$",
                    x = x[, 1L]
                ),
                allAreMatchingRegex(
                    pattern = "^ENS.*G[0-9]{11}\\.[0-9]+$",
                    x = x[, 2L]
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
                    pattern = "^FBtr[0-9]{7}$",
                    x = x[, 1L]
                ),
                allAreMatchingRegex(
                    pattern = "^FBgn[0-9]{7}$",
                    x = x[, 2L]
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
                    pattern = "^ENS.*T[0-9]{11}\\.[0-9]+(_PAR_Y)?$",
                    x = x[, 1L]
                ),
                allAreMatchingRegex(
                    pattern = "^ENS.*G[0-9]{11}\\.[0-9]+(_PAR_Y)?$",
                    x = x[, 2L]
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
            x[, 2L] <- gsub(pattern = "^gene=", replacement = "", x = x[, 2L])
            assert(
                allAreMatchingRegex(
                    pattern = "^[A-Za-z0-9_]+\\.[a-z0-9]+\\.[0-9]+$",
                    x = x[, 1L]
                ),
                allAreMatchingRegex(
                    pattern = "^WBGene[0-9]{8}$",
                    x = x[, 2L]
                )
            )
        },
        stop(sprintf("Unsupported FASTA: '%s'.", basename(file)))
    )
    out <- Tx2Gene(x)
    meta <- list(
        "call" = tryCatch(
            expr = standardizeCall(),
            error = function(e) NULL
        ),
        "date" = Sys.Date(),
        "file" = file,
        "packageVersion" = .pkgVersion,
        "provider" = provider
    )
    metadata(out) <- meta
    out
}
