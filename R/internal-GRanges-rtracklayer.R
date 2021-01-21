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
    genes <- object
    if (level == "transcripts") {
        transcripts <- object
    }
    rm(object)
    ## Genes -------------------------------------------------------------------
    ## These annotations will be included at transcript level (see below).
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
        out <- genes
    }
    ## Transcripts -------------------------------------------------------------
    if (level == "transcripts") {
        transcripts <-
            transcripts[!is.na(sanitizeNA(mcols(transcripts)[["transcript_id"]]))]
        assert(hasLength(transcripts))
        if (source == "Ensembl" && type == "GFF3") {
            transcripts <- .makeTranscriptsFromEnsemblGFF3(transcripts)
        } else if (source == "Ensembl" && type == "GTF") {
            transcripts <- .makeTranscriptsFromEnsemblGTF(transcripts)
        } else if (source == "FlyBase" && type == "GTF") {
            transcripts <- .makeTranscriptsFromFlyBaseGTF(transcripts)
        } else if (source == "GENCODE" && type == "GFF3") {
            transcripts <- .makeTranscriptsFromGencodeGFF3(transcripts)
        } else if (source == "GENCODE" && type == "GTF") {
            transcripts <- .makeTranscriptsFromGencodeGTF(transcripts)
        } else if (source == "RefSeq" && type == "GFF3") {
            transcripts <- .makeTranscriptsFromRefSeqGFF3(transcripts)
        } else if (source == "RefSeq" && type == "GTF") {
            transcripts <- .makeTranscriptsFromRefSeqGTF(transcripts)
        } else if (source == "WormBase" && type == "GTF") {
            transcripts <- .makeTranscriptsFromWormBaseGTF(transcripts)
        } else {
            ## nocov start
            stop(
                "Failed to make transcript-level GRanges.\n",
                "Unsupported GFF file format."
            )
            ## nocov end
        }
        ## Remove GFF-specific parent columns, etc.
        if (type == "GFF3") {
            transcripts <- .minimizeGFF3(transcripts)
        }
        ## Set names and stash metadata.
        names(transcripts) <- mcols(transcripts)[["transcript_id"]]
        metadata(transcripts)[["level"]] <- "transcripts"
        ## Skip gene-level metadata merge for GRanges that have been split
        ## into GRangesList.
        if (
            is(genes, "GRanges") &&
            ## This step is necessary for RefSeq GFF3.
            !anyDuplicated(mcols(genes)[["gene_id"]])
        ) {
            ## By default, merge the gene-level annotations into the
            ## transcript-level ones, for objects that have ranges 1:1 with the
            ## identifiers.
            out <- .mergeGenesIntoTranscripts(transcripts, genes)
        } else {
            cli_alert_warning("Skipping gene metadata merge.")
            out <- transcripts
        }
    }
}
