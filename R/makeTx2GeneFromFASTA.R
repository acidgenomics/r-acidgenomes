## FIXME THIS NEEDS TO SUPPORT IGNORE VERSION ARGUMENT.
## FIXME DETECT THE SOURCE HERE AUTOMATICALLY.



#' Make a Tx2Gene object from transcriptome FASTA
#'
#' @export
#' @note RefSeq transcript FASTA
#'   (e.g. "GCF_000001405.39_GRCh38.p13_rna.fna.gz") doesn't contain gene
#'   identifiers, and is not supported.
#' @note Updated 2021-01-21.
#'
#' @inheritParams AcidRoxygen::params
#' @param source `character(1)`.
#'   FASTA file source:
#'
#'   - `"ensembl"`: Ensembl.
#'   - `"gencode"`: GENCODE.
#'   - `"flybase"`: FlyBase.
#'   - `"wormbase"`: WormBase.
#'
#'   Assuming Ensembl transcriptome (i.e. cDNA) input by default.
#'
#' @return `Tx2Gene`.
#'
#' @examples
#' ## Ensembl ====
#' file <- pasteURL(
#'     "ftp.ensembl.org",
#'     "pub",
#'     "release-98",
#'     "fasta",
#'     "homo_sapiens",
#'     "cdna",
#'     "Homo_sapiens.GRCh38.cdna.all.fa.gz",
#'     protocol = "ftp"
#' )
#' ## > makeTx2GeneFromFASTA(file, source = "ensembl")
#'
#' ## GENCODE ====
#' file <- pasteURL(
#'     "ftp.ebi.ac.uk",
#'     "pub",
#'     "databases",
#'     "gencode",
#'     "Gencode_human",
#'     "release_32",
#'     "gencode.v32.transcripts.fa.gz",
#'     protocol = "ftp"
#' )
#' ## > makeTx2GeneFromFASTA(file, source = "gencode")
#'
#' ## FlyBase ====
#' file <- pasteURL(
#'     "ftp.flybase.net",
#'     "releases",
#'     "FB2019_05",
#'     "dmel_r6.30",
#'     "fasta",
#'     "dmel-all-transcript-r6.30.fasta.gz",
#'     protocol = "ftp"
#' )
#' ## > makeTx2GeneFromFASTA(file, source = "flybase")
#'
#' ## WormBase ====
#' file <- pasteURL(
#'     "ftp.wormbase.org",
#'     "pub",
#'     "wormbase",
#'     "releases",
#'     "WS272",
#'     "species",
#'     "c_elegans",
#'     "PRJNA13758",
#'     "c_elegans.PRJNA13758.WS272.mRNA_transcripts.fa.gz",
#'     protocol = "ftp"
#' )
#' ## > makeTx2GeneFromFASTA(file, source = "wormbase")
makeTx2GeneFromFASTA <- function(
    file,
    source = c(
        "ensembl",
        "gencode",
        "flybase",
        "wormbase"
    )
) {
    file <- .cacheIt(file)
    x <- import(file, format = "lines")
    source <- match.arg(source)
    x <- grep(pattern = "^>", x = x, value = TRUE)
    if (!hasLength(x)) {
        stop("FASTA file does not contain '>' annotations.")
    }
    x <- substr(x, start = 2L, stop = nchar(x))
    switch(
        EXPR = source,
        "ensembl" = {
            x <- strsplit(x = x, split = " ", fixed = TRUE)
            x <- lapply(
                X = x,
                FUN = function(x) {
                    x[c(1L, 4L)]
                }
            )
            x <- do.call(what = rbind, args = x)
            x[, 2L] <- gsub(pattern = "^gene:", replacement = "", x = x[, 2L])
        },
        "gencode" = {
            x <- strsplit(x = x, split = "|", fixed = TRUE)
            x <- lapply(
                X = x,
                FUN = function(x) {
                    x[c(1L, 2L)]
                }
            )
            x <- do.call(what = rbind, args = x)
        },
        "flybase" = {
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
        },
        "wormbase" = {
            x <- strsplit(x = x, split = " ", fixed = TRUE)
            x <- lapply(
                X = x,
                FUN = function(x) {
                    x[c(1L, 2L)]
                }
            )
            x <- do.call(what = rbind, args = x)
            x[, 2L] <- gsub(pattern = "^gene=", replacement = "", x = x[, 2L])
        }
    )
    x <- unique(x)
    Tx2Gene(x)
}
