.makeGRangesFromRtracklayer <- function(
    object,
    level = c("genes", "transcripts"),
    source,
    type
) {
    assert(is(object, "GRanges"))
    level <- match.arg(level)
    source <- match.arg(
        arg = source,
        choices = c("Ensembl", "FlyBase", "GENCODE", "RefSeq", "WormBase")
    )
    type <- match.arg(arg = type, choices = c("GFF3", "GTF"))
    ## Standardize -------------------------------------------------------------
    ## Standardize FlyBase, GENCODE, and RefSeq files to follow expected
    ## Ensembl-like naming conventions.
    object <- switch(
        EXPR = source,
        "FlyBase" = .standardizeFlyBaseToEnsembl(object),
        "GENCODE" = .standardizeGencodeToEnsembl(object),
        "RefSeq"  = .standardizeRefSeqToEnsembl(object),
        object
    )
    assert(
        isSubset(
            x = c("gene_id", "transcript_id"),
            y = colnames(mcols(object))
        ),
        ## `gene_type` needs to be renamed to `gene_biotype`, if defined.
        areDisjointSets(
            x = c("gene", "gene_type"),
            y = colnames(mcols(object))
        )
    )
    ## Genes -------------------------------------------------------------------
    ## These annotations will be included at transcript level (see below).
    genes <- object
    genes <- genes[!is.na(sanitizeNA(mcols(genes)[["gene_id"]]))]
    genes <- genes[is.na(sanitizeNA(mcols(genes)[["transcript_id"]]))]
    assert(hasLength(genes))
    what <- switch(
        EXPR = type,
        "GFF3" = {
            switch(
                EXPR = source,
                "Ensembl" = .makeGenesFromEnsemblGFF3,
                "GENCODE" = .makeGenesFromGencodeGFF3,
                "RefSeq" = .makeGenesFromRefSeqGFF3,
                NULL
            )
        },
        "GTF" = {
            switch(
                EXPR = source,
                "Ensembl" = .makeGenesFromEnsemblGTF,
                "FlyBase" = .makeGenesFromFlyBaseGTF,
                "GENCODE" = .makeGenesFromGencodeGTF,
                "RefSeq" = .makeGenesFromRefSeqGTF,
                "WormBase" = .makeGenesFromWormBaseGTF,
                NULL
            )
        }
    )
    if (!is.function(what)) {
        stop(sprintf("Unsupported genome file: %s %s.", source, type))
    }
    genes <- do.call(what = what, args = list(object = genes))
    ## Remove GFF-specific parent columns, etc.
    if (type == "GFF3") {
        genes <- .minimizeGFF3(genes)
    }
    mcols(genes) <- removeNA(mcols(genes))
    if (level == "genes") {
        return(genes)
    }
    ## Transcripts -------------------------------------------------------------
    tx <- object
    colnames(mcols(tx)) <-
        gsub(
            pattern = "^transcript_",
            replacement = "tx_",
            x = colnames(mcols(tx))
        )
    tx <- tx[!is.na(sanitizeNA(mcols(tx)[["tx_id"]]))]
    assert(hasLength(tx))
    what <- switch(
        EXPR = type,
        "GFF3" = switch(
            EXPR = source,
            "Ensembl" = .makeTranscriptsFromEnsemblGFF3,
            "GENCODE" = .makeTranscriptsFromGencodeGFF3,
            "RefSeq" = .makeTranscriptsFromRefSeqGFF3
        ),
        "GTF" = switch(
            EXPR = source,
            "Ensembl" = .makeTranscriptsFromEnsemblGTF,
            "FlyBase" = .makeTranscriptsFromFlyBaseGTF,
            "GENCODE" = .makeTranscriptsFromGencodeGTF,
            "RefSeq" = .makeTranscriptsFromRefSeqGTF,
            "WormBase" = .makeTranscriptsFromWormBaseGTF
        )
    )
    if (!is.function(what)) {
        stop(sprintf("Unsupported genome file: %s %s.", source, type))
    }
    tx <- do.call(what = what, args = list(object = tx))
    ## Remove GFF-specific parent columns, etc.
    if (type == "GFF3") {
        tx <- .minimizeGFF3(tx)
    }
    tx <- .mergeGenesIntoTranscripts(tx, genes)
    tx
}
