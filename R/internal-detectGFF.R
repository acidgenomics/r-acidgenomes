#' Detect the GFF source information
#'
#' @details
#' Assuming we've already cached the URL using BiocFileCache here.
#' This step will load into GRanges via rtracklayer.
#'
#' @note Updated 2021-01-09.
#' @noRd
.detectGFF <- function(file) {
    assert(isAFile(file))
    alert("Detecting annotation source.")
    gr <- import(file)
    assert(is(gr, "GRanges"))
    source <- .detectGFFSource(gr)
    type <- .detectGFFType(gr)
    alertInfo(sprintf("%s %s detected.", source, type))
    out <- c(
        "source" = source,
        "type" = type
    )
    out
}



#' Report the source of the gene annotations
#'
#' @note Updated 2021-01-10.
#' @noRd
.detectGFFSource <- function(object) {
    assert(is(object, "GRanges"))
    mcols <- mcols(object)
    source <- mcols[["source"]]
    if (
        ## UCSC (e.g. hg38_knownGene)
        any(grepl(pattern = "_knownGene$", x = source, ignore.case = FALSE))
    ) {
        ## nocov start
        stop(paste0(
            "UCSC genomes are intentionally not supported.\n",
            "Use a pre-built TxDb package instead ",
            "(e.g. 'TxDb.Hsapiens.UCSC.hg38.knownGene')."
        ))
        ## nocov end
    } else if (
        ## Check for GENCODE prior to Ensembl.
        any(source == "ENSEMBL") &&
        any(source == "HAVANA") &&
        "gene_type" %in% colnames(mcols)
    ) {
        out <- "GENCODE"
    } else if (
        any(grepl(pattern = "FlyBase", x = source, ignore.case = FALSE))
    ) {
        out <- "FlyBase"
    } else if (
        any(grepl(pattern = "WormBase", x = source, ignore.case = FALSE))
    ) {
        out <- "WormBase"
    } else if (
        any(grepl(pattern = "RefSeq", x = source, ignore.case = FALSE))
    ) {
        out <- "RefSeq"
    } else if (
        any(grepl(
            pattern = "ensembl|havana",
            x = source,
            ignore.case = FALSE
        ))
    ) {
        out <- "Ensembl"
    } else {
        ## nocov start
        stop(sprintf(
            fmt = paste(
                "Failed to detect valid GFF/GTF source.",
                "Supported: %s",
                sep = "\n"
            ),
            toString(c("Ensembl", "FlyBase", "GENCODE", "RefSeq", "WormBase"))
        ))
        ## nocov end
    }
    out
}



#' Determine if GFF or GTF
#'
#' @note Updated 2021-01-10.
#' @noRd
.detectGFFType <- function(object) {
    assert(is(object, "GRanges"))
    if (any(c("ID", "Name", "Parent") %in% colnames(mcols(object)))) {
        out <- "GFF3"
    } else {
        out <- "GTF"
    }
    out
}
