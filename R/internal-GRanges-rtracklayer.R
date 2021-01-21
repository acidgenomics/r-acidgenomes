.makeGRangesFromRtracklayer <- function(
    object,
    level = c("genes", "transcripts"),
    source
) {
    assert(is(object, "GRanges"))
    level <- match.arg(level)
    source <- match.arg(
        arg = source,
        choices = c("Ensembl", "FlyBase", "GENCODE", "RefSeq", "WormBase")
    )
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
    mcolnames <- colnames(mcols(object))
    assert(
        isSubset(x = c("gene_id", "transcript_id"), y = mcolnames),
        ## `gene_type` needs to be renamed to `gene_biotype`, if defined.
        areDisjointSets(x = c("gene", "gene_type"), y = mcolnames)
    )
    rm(mcolnames)
    genes <- object
    if (level == "transcripts") {
        transcripts <- object
    }
    rm(object)
    ## Genes -------------------------------------------------------------------
    ## `makeGRangesFromGFF()` attempts to always returns gene-level metadata,
    ## even when transcripts are requested. We'll merge this object into the
    ## transcript-level GRanges below, when possible.
    genes <- genes[!is.na(mcols(genes)[["gene_id"]])]
    genes <- genes[is.na(mcols(genes)[["transcript_id"]])]
    assert(hasLength(genes))
    if (source == "Ensembl" && type == "GFF3") {
        genes <- .makeGenesFromEnsemblGFF3(genes)
    } else if (source == "Ensembl" && type == "GTF") {
        genes <- .makeGenesFromEnsemblGTF(genes)
    } else if (source == "FlyBase" && type == "GTF") {
        genes <- .makeGenesFromFlyBaseGTF(genes)
    } else if (source == "GENCODE" && type == "GFF3") {
        genes <- .makeGenesFromGencodeGFF3(genes)
    } else if (source == "GENCODE" && type == "GTF") {
        genes <- .makeGenesFromGencodeGTF(genes)
    } else if (source == "RefSeq" && type == "GFF3") {
        genes <- .makeGenesFromRefSeqGFF3(genes)
    } else if (source == "RefSeq" && type == "GTF") {
        genes <- .makeGenesFromRefSeqGTF(genes)
    } else if (source == "WormBase" && type == "GTF") {
        genes <- .makeGenesFromWormBaseGTF(genes)
    } else {
        ## nocov start
        stop(
            "Failed to make gene-level GRanges.\n",
            "Unsupported GFF source file."
        )
        ## nocov end
    }
    ## Remove GFF-specific parent columns, etc.
    if (type == "GFF3") {
        genes <- .minimizeGFF3(genes)
    }
    ## Set names and stash metadata.
    names(genes) <- mcols(genes)[["gene_id"]]
    metadata(genes)[["level"]] <- "genes"
    if (level == "genes") {
        out <- genes
    }
    ## Transcripts -------------------------------------------------------------
    if (level == "transcripts") {
        transcripts <-
            transcripts[!is.na(mcols(transcripts)[["transcript_id"]])]
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
