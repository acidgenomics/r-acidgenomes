#' Generate empty genomic ranges
#'
#' Utility function that provides support for creating internal `GRanges` for
#' transgene and FASTA spike-in sequences.
#'
#' @export
#' @note Updated 2019-08-21.
#'
#' @param names `character`.
#' Gene or transcript names.
#'
#' @param seqname `character(1)`.
#' Name of the alternative chromosome to be defined in `seqnames` where these
#' ranges will be grouped. Defaults to `"unknown"` but `"transgene"`
#' (transgenes) and `"spike"` (spike-ins) are also supported.
#'
#' @param mcolnames `character` or `NULL`.
#' Metadata column names to be defined in the `mcols` of the `GRanges` return.
#' Normally this does not need to be defined; useful when combining with another
#' `GRanges` that contains metadata.
#'
#' @return `GRanges`.
#'
#' @seealso `help("seqinfo", "GenomeInfoDb")`.
#'
#' @examples
#' ## Unknown/dead genes.
#' emptyRanges("ENSG00000000000", seqname = "unknown")
#'
#' ## Transgenes.
#' emptyRanges(c("EGFP", "TDTOMATO", "GAL4"), seqname = "transgene")
#'
#' ## Spike-ins.
#' emptyRanges("ERCC", seqname = "spike")
emptyRanges <-
    function(
        names,
        seqname = c("unknown", "transgene", "spike"),
        mcolnames = NULL
    ) {
        assert(
            isCharacter(names),
            isAny(mcolnames, c("character", "NULL"))
        )
        seqname <- match.arg(seqname)
        gr <- GRanges(
            seqnames = seqname,
            ranges = IRanges(
                start = (seq_len(length(names)) - 1L) * 100L + 1L,
                width = 100L
            )
        )
        names(gr) <- names
        ## Create the required empty metadata columns.
        if (!hasLength(mcolnames)) {
            ncol <- 0L
        } else {
            ncol <- length(mcolnames)
        }
        mcols <- matrix(
            nrow = length(names),
            ncol = ncol,
            dimnames = list(names, mcolnames)
        )
        mcols <- as(mcols, "DFrame")
        mcols(gr) <- mcols
        gr
    }
