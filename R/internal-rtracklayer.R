#' Determine if input is GFF (GFF3) or GTF (GFFv2)
#'
#' @note Updated 2021-01-25.
#' @noRd
.rtracklayerFormat <- function(object) {
    assert(is(object, "GRanges"))
    if (any(c("ID", "Name", "Parent") %in% colnames(mcols(object)))) {
        x <- "GFF"
    } else {
        x <- "GTF"
    }
    x
}



#' Detect the provider (i.e. source) of the genome annotations
#'
#' @note Updated 2021-01-25.
#' @noRd
.rtracklayerProvider <- function(object) {
    assert(is(object, "GRanges"))
    mcols <- mcols(object)
    source <- mcols[["source"]]
    if (
        ## UCSC (e.g. hg38_knownGene).
        any(grepl(
            pattern = "_(ensGene|knownGene|ncbiRefSeq|refGene)$",
            x = source,
            ignore.case = FALSE
        ))
    ) {
        x <- "UCSC"
    } else if (
        ## Check for GENCODE prior to Ensembl.
        any(source == "ENSEMBL") &&
        any(source == "HAVANA") &&
        "gene_type" %in% colnames(mcols)
    ) {
        x <- "GENCODE"
    } else if (
        any(grepl(pattern = "FlyBase", x = source, ignore.case = FALSE))
    ) {
        x <- "FlyBase"
    } else if (
        any(grepl(pattern = "WormBase", x = source, ignore.case = FALSE))
    ) {
        x <- "WormBase"
    } else if (
        any(grepl(pattern = "RefSeq", x = source, ignore.case = FALSE))
    ) {
        x <- "RefSeq"
    } else if (
        any(grepl(
            pattern = "ensembl|havana",
            x = source,
            ignore.case = FALSE
        ))
    ) {
        x <- "Ensembl"
    } else {
        ## nocov start
        stop(sprintf(
            fmt = paste(
                "Failed to detect valid GFF/GTF source.",
                "Supported: %s",
                sep = "\n"
            ),
        ))
        ## nocov end
    }
    x
}



## FIXME Is it possible to return this without the gene-level merge step?

## Updated 2021-01-24.
.makeGRangesFromRtracklayer <- function(
    file,
    level = c("genes", "transcripts")
) {
    object <- import(file = .cacheIt(file))
    assert(is(object, "GRanges"))
    level <- match.arg(level)
    format <- match.arg(
        arg = .grangesFormat(object),
        choices = c("GFF", "GTF")
    )
    provider <- match.arg(
        arg = .grangesProvider(object),
        choices = c(
            "Ensembl",
            "FlyBase",
            "GENCODE",
            "RefSeq",
            "UCSC",
            "WormBase"
        )
    )
    ## Standardize -------------------------------------------------------------
    ## Standardize FlyBase, GENCODE, and RefSeq files to follow expected
    ## Ensembl-like naming conventions.
    object <- switch(
        EXPR = provider,
        "FlyBase" = .standardizeFlyBaseToEnsembl(object),
        "GENCODE" = .standardizeGencodeToEnsembl(object),
        "RefSeq"  = .standardizeRefSeqToEnsembl(object),
        object
    )
    colnames(mcols(object)) <-
        gsub(
            pattern = "^transcript_",
            replacement = "tx_",
            x = colnames(mcols(object))
        )
    assert(
        isSubset(
            x = c("gene_id", "tx_id"),
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
    genes <- genes[is.na(sanitizeNA(mcols(genes)[["tx_id"]]))]
    assert(hasLength(genes))
    what <- switch(
        EXPR = format,
        "GFF" = {
            switch(
                EXPR = provider,
                "Ensembl" = .makeGenesFromEnsemblGFF,
                "GENCODE" = .makeGenesFromGencodeGFF,
                "RefSeq" = .makeGenesFromRefSeqGFF,
                NULL
            )
        },
        "GTF" = {
            switch(
                EXPR = provider,
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
        stop(sprintf("Unsupported genome file: %s %s.", provider, type))
    }
    genes <- do.call(what = what, args = list("object" = genes))
    mcols(genes) <- removeNA(mcols(genes))
    geneCol <- .matchGRangesNamesColumn(genes)
    assert(hasNoDuplicates(mcols(genes)[[geneCol]]))
    if (level == "genes") {
        return(genes)
    }
    ## Transcripts -------------------------------------------------------------
    tx <- object
    tx <- tx[!is.na(sanitizeNA(mcols(tx)[["tx_id"]]))]
    assert(hasLength(tx))
    what <- switch(
        EXPR = format,
        "GFF" = switch(
            EXPR = provider,
            "Ensembl" = .makeTranscriptsFromEnsemblGFF,
            "GENCODE" = .makeTranscriptsFromGencodeGFF,
            "RefSeq" = .makeTranscriptsFromRefSeqGFF
        ),
        "GTF" = switch(
            EXPR = provider,
            "Ensembl" = .makeTranscriptsFromEnsemblGTF,
            "FlyBase" = .makeTranscriptsFromFlyBaseGTF,
            "GENCODE" = .makeTranscriptsFromGencodeGTF,
            "RefSeq" = .makeTranscriptsFromRefSeqGTF,
            "WormBase" = .makeTranscriptsFromWormBaseGTF
        )
    )
    if (!is.function(what)) {
        stop(sprintf("Unsupported genome file: %s %s.", provider, type))
    }
    tx <- do.call(what = what, args = list(object = tx))
    mcols(tx) <- removeNA(mcols(tx))
    tx <- .mergeGenesIntoTranscripts(tx, genes)
    tx
}



## Ensembl =====================================================================

## GTF:
## >  [1] "source"                   "type"
## >  [3] "score"                    "phase"
## >  [5] "gene_id"                  "gene_version"
## >  [7] "gene_name"                "gene_source"
## >  [9] "gene_biotype"             "transcript_id"
## > [11] "transcript_version"       "transcript_name"
## > [13] "transcript_source"        "transcript_biotype"
## > [15] "tag"                      "transcript_support_level"
## > [17] "exon_number"              "exon_id"
## > [19] "exon_version"             "protein_id"
## > [21] "protein_version"          "ccds_id"
##
## GFF:
## >  [1] "source"                   "type"
## >  [3] "score"                    "phase"
## >  [5] "ID"                       "Alias"
## >  [7] "external_name"            "logic_name"
## >  [9] "Name"                     "biotype"
## > [11] "description"              "gene_id"
## > [13] "version"                  "Parent"
## > [15] "tag"                      "transcript_id"
## > [17] "transcript_support_level" "constitutive"
## > [19] "ensembl_end_phase"        "ensembl_phase"
## > [21] "exon_id"                  "rank"
## > [23] "protein_id"               "ccdsid"



## Note that call upstream in `.makeGenesFromGFF()` will prepare the rows
## properly already, by filtering aganist `gene_id` and `transcript_id`.
.makeGenesFromEnsemblGFF <- function(object) {
    assert(is(object, "GRanges"))
    ## Assign `gene_name` from `Name` column.
    assert(
        isSubset("Name", colnames(mcols(object))),
        areDisjointSets("gene_name", colnames(mcols(object)))
    )
    mcols(object)[["gene_name"]] <- mcols(object)[["Name"]]
    ## Assign `gene_biotype` from `biotype` column.
    assert(
        isSubset("biotype", colnames(mcols(object))),
        areDisjointSets("gene_biotype", colnames(mcols(object)))
    )
    mcols(object)[["gene_biotype"]] <- mcols(object)[["biotype"]]
    object
}



.makeGenesFromEnsemblGTF <- function(object) {
    assert(
        is(object, "GRanges"),
        isSubset(
            x = c("gene_id", "type"),
            y = colnames(mcols(object))
        )
    )
    object <- object[mcols(object)[["type"]] == "gene"]
    object
}



.makeTranscriptsFromEnsemblGFF <- function(object) {
    assert(is(object, "GRanges"))
    ## Assign `transcript_name` from `Name` column.
    assert(
        isSubset("Name", colnames(mcols(object))),
        areDisjointSets("transcript_name", colnames(mcols(object)))
    )
    mcols(object)[["transcript_name"]] <- mcols(object)[["Name"]]
    ## Assign `transcript_biotype` from `biotype` column.
    assert(
        isSubset("biotype", colnames(mcols(object))),
        areDisjointSets("transcript_biotype", colnames(mcols(object)))
    )
    mcols(object)[["transcript_biotype"]] <- mcols(object)[["biotype"]]
    ## Assign `gene_id` from `Parent` column.
    assert(
        isSubset("Parent", colnames(mcols(object))),
        all(grepl("^gene:", mcols(object)[["Parent"]]))
    )
    mcols(object)[["gene_id"]] <- as.character(mcols(object)[["Parent"]])
    mcols(object)[["gene_id"]] <- gsub(
        pattern = "^gene:",
        replacement = "",
        x = mcols(object)[["gene_id"]]
    )
    object
}



.makeTranscriptsFromEnsemblGTF <- function(object) {
    assert(
        is(object, "GRanges"),
        isSubset(
            x = c("transcript_id", "type"),
            y = colnames(mcols(object))
        )
    )
    object <- object[mcols(object)[["type"]] == "transcript"]
    object
}



## FlyBase =====================================================================

## GTF:
## > [1] "source"            "type"              "score"
## > [4] "phase"             "gene_id"           "gene_symbol"
## > [7] "transcript_id"     "transcript_symbol" "#"

## Compatible with Ensembl importer after we run `.standardizeFlyBaseGFF()`,
## which is called in `.makeGenesFromGFF()`.
.makeGenesFromFlyBaseGTF <- function(object) {
    assert(is(object, "GRanges"))
    object <- .makeGenesFromEnsemblGTF(object)
    object
}



.makeTranscriptsFromFlyBaseGTF <- function(object) {
    assert(is(object, "GRanges"))
    ## Note that FlyBase uses non-standard transcript types.
    keep <- grepl(
        pattern = paste(c("^pseudogene$", "RNA$"), collapse = "|"),
        x = mcols(object)[["type"]],
        ignore.case = TRUE
    )
    object <- object[keep]
    object
}



.standardizeFlyBaseToEnsembl <- function(object) {
    assert(is(object, "GRanges"))
    mcolnames <- colnames(mcols(object))
    ## Match Ensembl spec by renaming `*_symbol` to `*_name`.
    mcolnames <- sub(
        pattern = "^gene_symbol$",
        replacement = "gene_name",
        x = mcolnames
    )
    mcolnames <- sub(
        pattern = "^transcript_symbol$",
        replacement = "transcript_name",
        x = mcolnames
    )
    colnames(mcols(object)) <- mcolnames
    object
}



## GENCODE =====================================================================

## Uses `gene_type` instead of `gene_biotype`.
## Note that `gene_id` and `gene_name` are nicely defined, so don't use `Name`.
## Consider removing gene and transcript versions automatically.
##
## GTF:
## >  [1] "source"                   "type"
## >  [3] "score"                    "phase"
## >  [5] "gene_id"                  "gene_type"
## >  [7] "gene_name"                "level"
## >  [9] "havana_gene"              "transcript_id"
## > [11] "transcript_type"          "transcript_name"
## > [13] "transcript_support_level" "tag"
## > [15] "havana_transcript"        "exon_number"
## > [17] "exon_id"                  "ont"
## > [19] "protein_id"               "ccdsid"
##
## GFF:
## >  [1] "source"                   "type"
## >  [3] "score"                    "phase"
## >  [5] "ID"                       "gene_id"
## >  [7] "gene_type"                "gene_name"
## >  [9] "level"                    "havana_gene"
## > [11] "Parent"                   "transcript_id"
## > [13] "transcript_type"          "transcript_name"
## > [15] "transcript_support_level" "tag"
## > [17] "havana_transcript"        "exon_number"
## > [19] "exon_id"                  "ont"
## > [21] "protein_id"               "ccdsid"



#' Detect PAR duplicates
#'
#' Match Ensembl spec by removing the duplicate PAR Y chromosome annotations.
#'
#' @note Updated 2021-01-24.
#' @noRd
.detectPARDupes <- function(object, idCol) {
    assert(is(object, "GRanges"))
    idCol <- match.arg(
        arg = idCol,
        choices = c("ID", "gene_id", "transcript_id")
    )
    dupes <- grep(pattern = "_PAR_Y$", x = mcols(object)[[idCol]], value = TRUE)
    if (hasLength(dupes)) {
        alertWarning(sprintf(
            "%d pseudoautosomal region (PAR) Y chromosome %s: {.var %s}.",
            length(dupes),
            ngettext(
                n = length(dupes),
                msg1 = "duplicate",
                msg2 = "duplicates"
            ),
            toString(dupes, width = 100L)
        ))
    }
    invisible(object)
}



.makeGenesFromGencodeGFF <- function(object) {
    assert(
        is(object, "GRanges"),
        isSubset(
            x = c("ID", "gene_biotype", "gene_id", "gene_name", "type"),
            y = colnames(mcols(object))
        )
    )
    object <- object[mcols(object)[["type"]] == "gene"]
    .detectPARDupes(object, idCol = "ID")
    object
}



.makeGenesFromGencodeGTF <- function(object) {
    object <- .makeGenesFromEnsemblGTF(object)
    .detectPARDupes(object, idCol = "gene_id")
    object
}



.makeTranscriptsFromGencodeGFF <- function(object) {
    assert(
        is(object, "GRanges"),
        isSubset(
            x = c(
                "gene_biotype", "gene_id",
                "transcript_biotype", "transcript_id"
            ),
            y = colnames(mcols(object))
        )
    )
    object <- object[mcols(object)[["type"]] == "transcript"]
    .detectPARDupes(object, idCol = "ID")
    object
}



.makeTranscriptsFromGencodeGTF <- function(object) {
    object <- .makeTranscriptsFromEnsemblGTF(object)
    .detectPARDupes(object, idCol = "transcript_id")
    object
}



.standardizeGencodeToEnsembl <- function(object) {
    assert(is(object, "GRanges"))
    mcolnames <- colnames(mcols(object))
    ## Match Ensembl spec, which uses `*_biotype` instead of `*_type`.
    mcolnames <- sub(
        pattern = "^gene_type$",
        replacement = "gene_biotype",
        x = mcolnames
    )
    mcolnames <- sub(
        pattern = "^transcript_type$",
        replacement = "transcript_biotype",
        x = mcolnames
    )
    colnames(mcols(object)) <- mcolnames
    object
}



## RefSeq ======================================================================

## Note that GTF contains both "gene_id" and "gene" columns, whereas GFF3 format
## only contains "gene" column.
##
## GTF:
## >  [1] "source"             "type"               "score"
## >  [4] "phase"              "gene_id"            "db_xref"
## >  [7] "description"        "gbkey"              "gene"
## > [10] "gene_biotype"       "pseudo"             "transcript_id"
## > [13] "product"            "exon_number"        "gene_synonym"
## > [16] "model_evidence"     "protein_id"         "exception"
## > [19] "inference"          "note"               "anticodon"
## > [22] "partial"            "transl_except"      "standard_name"
## > [25] "ribosomal_slippage" "codons"             "transl_table"
##
## GFF:
## >  [1] "source"                    "type"
## >  [3] "score"                     "phase"
## >  [5] "ID"                        "Dbxref"
## >  [7] "Name"                      "chromosome"
## >  [9] "gbkey"                     "genome"
## > [11] "mol_type"                  "description"
## > [13] "gene"                      "gene_biotype"
## > [15] "pseudo"                    "Parent"
## > [17] "product"                   "transcript_id"
## > [19] "gene_synonym"              "model_evidence"
## > [21] "protein_id"                "Note"
## > [23] "exception"                 "inference"
## > [25] "standard_name"             "experiment"
## > [27] "function"                  "regulatory_class"
## > [29] "feat_class"                "recombination_class"
## > [31] "rpt_type"                  "rpt_unit_seq"
## > [33] "anticodon"                 "partial"
## > [35] "start_range"               "end_range"
## > [37] "transl_except"             "mobile_element_type"
## > [39] "rpt_family"                "satellite"
## > [41] "bound_moiety"              "Target"
## > [43] "assembly_bases_aln"        "assembly_bases_seq"
## > [45] "bit_score"                 "blast_aligner"
## > [47] "blast_score"               "common_component"
## > [49] "e_value"                   "filter_score"
## > [51] "for_remapping"             "gap_count"
## > [53] "hsp_percent_coverage"      "matchable_bases"
## > [55] "matched_bases"             "num_ident"
## > [57] "num_mismatch"              "pct_coverage"
## > [59] "pct_coverage_hiqual"       "pct_identity_gap"
## > [61] "pct_identity_gapopen_only" "pct_identity_ungap"
## > [63] "rank"                      "weighted_identity"
## > [65] "lxr_locAcc_currStat_120"   "not_for_annotation"
## > [67] "consensus_splices"         "exon_identity"
## > [69] "identity"                  "idty"
## > [71] "matches"                   "product_coverage"
## > [73] "splices"                   "Gap"
## > [75] "merge_aligner"             "map"
## > [77] "part"                      "lxr_locAcc_currStat_35"
## > [79] "direction"                 "rpt_unit_range"
## > [81] "exon_number"               "number"
## > [83] "allele"                    "align_id"
## > [85] "batch_id"                  "crc32"
## > [87] "curated_alignment"         "promoted_rank"
## > [89] "qtaxid"                    "Is_circular"
## > [91] "country"                   "isolation-source"
## > [93] "note"                      "tissue-type"
## > [95] "codons"                    "transl_table"
##
## Types that map to `gene_id`:
## > [1] "enhancer"             "gene"
## > [2] "promoter"             "pseudogene"
## > [5] "recombination_region" "sequence_feature"
##
## Types that map to `transcript_id`:
## >  [1] "antisense_RNA"      "exon"
## >  [3] "guide_RNA"          "lnc_RNA"
## >  [5] "mRNA"               "primary_transcript"
## >  [7] "RNase_MRP_RNA"      "RNase_P_RNA"
## >  [9] "rRNA"               "scRNA"
## > [11] "snoRNA"             "snRNA"
## > [13] "telomerase_RNA"     "transcript"
## > [15] "vault_RNA"          "Y_RNA"



## Updated 2020-01-20.
.makeGenesFromRefSeqGFF <- function(object) {
    assert(
        is(object, "GRanges"),
        isSubset(
            x = c("Name", "gene_biotype", "gene_id", "type"),
            y = colnames(mcols(object))
        )
    )
    ## Only keep annotations that map to `Name` column.
    keep <- !is.na(mcols(object)[["Name"]])
    object <- object[keep]
    ## Drop rows that contain a `Parent` element.
    keep <- bapply(
        X = mcols(object)[["Parent"]],
        FUN = function(x) {
            identical(x, character(0L))
        }
    )
    object <- object[keep]
    ## Define `gene_name` from `gene_id`.
    mcols(object)[["gene_name"]] <- mcols(object)[["gene_id"]]
    object
}



## Updated 2020-01-20.
.makeGenesFromRefSeqGTF <- function(object) {
    assert(
        is(object, "GRanges"),
        isSubset(c("gene_biotype", "gene_id", "type"), colnames(mcols(object))),
        areDisjointSets("gene_name", colnames(mcols(object)))
    )
    ## Define `gene_name` from `gene_id`.
    mcols(object)[["gene_name"]] <- mcols(object)[["gene_id"]]
    object <- object[mcols(object)[["type"]] == "gene"]
    object
}



## Updated 2020-01-20.
.makeTranscriptsFromRefSeqGFF <- function(object) {
    assert(
        is(object, "GRanges"),
        isSubset("Name", colnames(mcols(object))),
        areDisjointSets("transcript_name", colnames(mcols(object)))
    )
    ## Only keep annotations that map to `Name` column.
    keep <- !is.na(mcols(object)[["Name"]])
    object <- object[keep]
    object
}



## Updated 2020-01-20.
.makeTranscriptsFromRefSeqGTF <- function(object) {
    assert(
        is(object, "GRanges"),
        isSubset(c("transcript_id", "type"), colnames(mcols(object))),
        areDisjointSets(
            x = c("gene_name", "transcript_biotype", "transcript_name"),
            y = colnames(mcols(object))
        )
    )
    ## Note that we're filtering by "exon" instead of "transcript" here.
    keep <- mcols(object)[["type"]] == "exon"
    assert(any(keep))
    n <- sum(keep, na.rm = TRUE)
    alertInfo(sprintf(
        "%d %s detected.",
        n, ngettext(n = n, msg1 = "exon", msg2 = "exons")
    ))
    object <- object[keep]
    ## Define `gene_name` from `gene_id`.
    mcols(object)[["gene_name"]] <- mcols(object)[["gene_id"]]
    ## Define `transcript_biotype` from `gene_biotype`.
    mcols(object)[["transcript_biotype"]] <- mcols(object)[["gene_biotype"]]
    ## Define `transcript_name` from `transcript_id`.
    mcols(object)[["transcript_name"]] <- mcols(object)[["transcript_id"]]
    object
}



## This step ensures that `gene_id` and `transcript_id` columns are defined.
## Updated 2020-01-20.
.standardizeRefSeqToEnsembl <- function(object) {
    assert(is(object, "GRanges"))
    mcols <- mcols(object)
    assert(
        isSubset(
            x = c("gene", "transcript_id"),
            y = colnames(mcols)
        ),
        areDisjointSets(
            x = c("gene_name", "transcript_name"),
            y = colnames(mcols)
        )
    )
    ## Ensure `gene_id` is defined.
    if (isTRUE(all(c("gene", "gene_id") %in% colnames(mcols)))) {
        ## Pick `gene_id` over `gene` column, if both are defined.
        ## This applies to GTF spec.
        keep <- setdiff(colnames(mcols), "gene")
        mcols <- mcols[keep]
    } else if ("gene" %in% colnames(mcols)) {
        ## Rename `gene` column to `gene_id`, matching Ensembl spec.
        ## This applies to GFF3 spec.
        colnames(mcols) <- sub("^gene$", "gene_id", colnames(mcols))
    }
    assert(isSubset(c("gene_id", "transcript_id"), colnames(mcols)))
    mcols(object) <- mcols
    object
}



## WormBase ====================================================================

## WormBase identifier fix may be needed. WormBase GTF currently imports
## somewhat malformed, and the gene identifiers require additional sanitization
## to return correctly. Look for rows containing "Gene:" and "Transcript:" in ID
## columns.
##
## GTF:
## >  [1] "source"             "type"               "score"
## >  [4] "phase"              "gene_id"            "gene_source"
## >  [7] "gene_biotype"       "transcript_id"      "transcript_source"
## > [10] "transcript_biotype" "exon_number"        "exon_id"
## > [13] "protein_id"



.makeGenesFromWormBaseGTF <- function(object) {
    assert(
        is(object, "GRanges"),
        isSubset(x = "gene_id", y = colnames(mcols(object))),
        areDisjointSets(x = "gene_name", y = colnames(mcols(object)))
    )
    ## Sanitize the `gene_id` column, which can contain some malformed entries
    ## prefixed with "Gene:". We want to keep these entries.
    mcols(object)[["gene_id"]] <- gsub(
        pattern = "^Gene:",
        replacement = "",
        x = mcols(object)[["gene_id"]]
    )
    ## Now safe to only keep rows that match "WBGene" in `gene_id`.
    ## Note that WormBase GTF currently contains some malformed "Transcript:"
    ## entries that we want to drop with this step.
    keep <- grepl(
        pattern = "^WBGene[[:digit:]]{8}$",
        x = mcols(object)[["gene_id"]]
    )
    object <- object[keep, , drop = FALSE]
    ## Process using Ensembl conventions.
    object <- .makeGenesFromEnsemblGTF(object)
    object
}



.makeTranscriptsFromWormBaseGTF <- function(object) {
    assert(
        is(object, "GRanges"),
        isSubset(
            x = c("gene_id", "transcript_id"),
            y = colnames(mcols(object))
        ),
        areDisjointSets(
            x = c("gene_name", "transcript_name"),
            y = colnames(mcols(object))
        )
    )
    ## Sanitize the `gene_id` column, which can contain some malformed entries
    ## prefixed with "Gene:". We want to keep these entries.
    mcols(object)[["gene_id"]] <- gsub(
        pattern = "^Gene:",
        replacement = "",
        x = mcols(object)[["gene_id"]]
    )
    ## Process using Ensembl conventions.
    object <- .makeTranscriptsFromEnsemblGTF(object)
    object
}
