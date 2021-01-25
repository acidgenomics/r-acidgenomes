## FIXME ## `gene_type` needs to be renamed to `gene_biotype`, if defined.
## FIXME NEED TO SANITIZE COLS FOR REFSEQ
## genes <- genes[!is.na(sanitizeNA(mcols(genes)[["gene_id"]]))]
## genes <- genes[is.na(sanitizeNA(mcols(genes)[["tx_id"]]))]
##
## FIXME ENSURE REFSEQ TRANSCRIPTS RETURN AS FLAT GRANGES OBJECT.
## FIXME TEST FLYBASE GFF AND WORMBASE GFF.
## FIXME INCLUDE GENEVERSION HERE IF POSSIBLE WHEN `IGNOREVERSION` = FALSE
## FIXME CURRENT RELEASE VERSION DOESNT SLOT ORGANISM HERE CORRECTLY.
##       RETHINK THAT FOR GFF.
## FIXME MAKE SURE FILE IS CORRECT URL, NOT TMPFILE BEFORE RELEASING.
## FIXME synonyms only works with Ensembl identifiers, consider making that more
##       clear in documentation.



## Updated 2021-01-24.
.makeGRangesFromRtracklayer <- function(
    file,
    level = c("genes", "transcripts"),
    ignoreVersion = TRUE
) {
    level <- match.arg(level)
    meta <- getGFFMetadata(file)
    meta[["level"]] <- level
    gr <- import(file = .cacheIt(file))
    assert(
        is(gr, "GRanges"),
        isString(meta[["format"]]),
        isString(meta[["provider"]])
    )
    format <- ifelse(
        test = grepl(pattern = "GTF", x = meta[["format"]]),
        yes = "GTF",
        no = "GFF"
    )
    provider <- match.arg(
        arg = meta[["provider"]],
        choices = c(
            "Ensembl",
            "FlyBase",
            "GENCODE",
            "RefSeq",
            "UCSC",
            "WormBase"
        )
    )
    funName <- paste0(
        ".",
        camelCase(
            object = paste("rtracklayer", level, "from", provider, format),
            strict = TRUE
        )
    )
    what <- .getFun(funName)
    gr <- do.call(what = what, args = list("object" = gr))
    metadata(gr) <- meta
    mcols(gr) <- removeNA(mcols(gr))
    if (identical(format, "GFF")) {
        ## Remove capitalized keys in mcols.
        keep <- !grepl(pattern = "^[A-Z]", x = colnames(mcols(gr)))
        mcols(gr) <- mcols(gr)[keep]
    }
    switch(
        EXPR = level,
        "genes" = {
            assert(allAreNotMatchingRegex(
                pattern = "^(transcript|tx)_",
                x = colnames(mcols(gr))
            ))
        },
        "transcripts" = {
            assert(isSubset("gene_id", colnames(mcols(gr))))
            colnames(mcols(gr)) <-
                gsub(
                    pattern = "^transcript_",
                    replacement = "tx_",
                    x = colnames(mcols(gr))
                )
        }
    )
    ## Remove any uninformative blacklisted columns.
    blacklistCols <- c(
        ## e.g. Ensembl GFF. Use "gene_biotype", "tx_biotype" instead.
        "biotype",
        ## e.g. Ensembl GFF: "havana_homo_sapiens". Not informative.
        "logic_name"
        ## FIXME Other values to consider:
        ## "biotype",
        ## "end_range",
        ## "exception",
        ## "gbkey",
        ## "partial",
        ## "pseudo",
        ## "source",
        ## "start_range",
        ## "transl_except"
    )
    keep <- !colnames(mcols(gr)) %in% blacklistCols
    mcols(gr) <- mcols(gr)[keep]
    idCol <- .matchGRangesNamesColumn(gr)
    assert(hasNoDuplicates(mcols(gr)[[idCol]]))
    names(gr) <- mcols(gr)[[idCol]]
    seqinfo <- .getSeqinfo(meta)
    if (is(seqinfo, "Seqinfo")) {
        seqinfo(gr) <- seqinfo[seqlevels(gr)]
    }
    gr <- .makeGRanges(
        object = gr,
        ignoreVersion = ignoreVersion
    )
    gr
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



## Updated 2021-01-25.
.rtracklayerGenesFromEnsemblGff <-
    function(object) {
        assert(
            is(object, "GRanges"),
            isSubset(
                x = c("Name", "biotype", "gene_id"),
                y = colnames(mcols(object))
            ),
            areDisjointSets(
                x = c("gene_biotype", "gene_name"),
                y = colnames(mcols(object))
            ),
            areDisjointSets("gene_biotype", colnames(mcols(object)))
        )
        keep <- !is.na(mcols(object)[["gene_id"]])
        object <- object[keep]
        mcols(object)[["gene_biotype"]] <- mcols(object)[["biotype"]]
        mcols(object)[["gene_name"]] <- mcols(object)[["Name"]]
        object
    }



## Updated 2021-01-25.
.rtracklayerGenesFromEnsemblGtf <-
    function(object) {
        assert(
            is(object, "GRanges"),
            isSubset(
                x = c("gene_id", "type"),
                y = colnames(mcols(object))
            )
        )
        keep <- mcols(object)[["type"]] == "gene"
        object <- object[keep]
        object
    }



## Updated 2021-01-25.
.rtracklayerTranscriptsFromEnsemblGff <-
    function(object) {
        assert(
            is(object, "GRanges"),
            isSubset(
                x = c("Name", "Parent", "biotype"),
                y = colnames(mcols(object))
            ),
            areDisjointSets(
                x = c("transcript_biotype", "transcript_name"),
                y = colnames(mcols(object))
            )
        )
        keep <- !is.na(mcols(object)[["transcript_id"]])
        object <- object[keep]
        assert(
            allAreMatchingRegex(
                pattern = "^gene:",
                x = as.character(mcols(object)[["Parent"]])
            )
        )
        mcols(object)[["transcript_biotype"]] <- mcols(object)[["biotype"]]
        mcols(object)[["transcript_name"]] <- mcols(object)[["Name"]]
        mcols(object)[["gene_id"]] <- gsub(
            pattern = "^gene:",
            replacement = "",
            x = as.character(mcols(object)[["Parent"]])
        )
        object
    }



## Updated 2021-01-25.
.rtracklayerTranscriptsFromEnsemblGtf <-
    function(object) {
        assert(
            is(object, "GRanges"),
            isSubset(
                x = c("transcript_id", "type"),
                y = colnames(mcols(object))
            )
        )
        keep <- mcols(object)[["type"]] == "transcript"
        object <- object[keep]
        object
    }



## FlyBase =====================================================================

## GTF:
## > [1] "source"            "type"              "score"
## > [4] "phase"             "gene_id"           "gene_symbol"
## > [7] "transcript_id"     "transcript_symbol" "#"



## Updated 2021-01-25.
.rtracklayerGenesFromFlyBaseGtf <-
    function(object) {
        object <- .standardizeFlyBaseToEnsembl(object)
        .rtracklayerGenesFromEnsemblGtf(object)
    }



## Updated 2021-01-25.
.rtracklayerTranscriptsFromFlyBaseGtf <-
    function(object) {
        object <- .standardizeFlyBaseToEnsembl(object)
        keep <- grepl(
            pattern = paste(c("^pseudogene$", "RNA$"), collapse = "|"),
            x = mcols(object)[["type"]],
            ignore.case = TRUE
        )
        object <- object[keep]
        object
    }



## Updated 2021-01-25.
.standardizeFlyBaseToEnsembl <-
    function(object) {
        assert(is(object, "GRanges"))
        colnames(mcols(object)) <- sub(
            pattern = "^gene_symbol$",
            replacement = "gene_name",
            x = colnames(mcols(object))
        )
        colnames(mcols(object)) <- sub(
            pattern = "^transcript_symbol$",
            replacement = "transcript_name",
            x = colnames(mcols(object))
        )
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



## FIXME CAN SET SEQINFO HERE.

## Updated 2021-01-25.
.makeGenesFromGencodeGff <-
    function(object) {
        assert(
            is(object, "GRanges"),
            isSubset(
                x = c("ID", "gene_biotype", "gene_id", "gene_name", "type"),
                y = colnames(mcols(object))
            )
        )
        keep <- mcols(object)[["type"]] == "gene"
        object <- object[keep]
        object
    }



## FIXME CAN SET SEQINFO HERE.

## Updated 2021-01-25.
.rtracklayerGenesFromGencodeGtf <-
    function(object) {
        .rtracklayerGenesFromEnsemblGtf(object)
    }



## FIXME CAN SET SEQINFO HERE.

## Updated 2021-01-25.
.rtracklayerTranscriptsFromGencodeGff <-
    function(object) {
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
        keep <- mcols(object)[["type"]] == "transcript"
        object <- object[keep]
        object
    }



## FIXME CAN SET SEQINFO HERE.

## Updated 2021-01-25.
.makeTranscriptsFromGencodeGtf <-
    function(object) {
        .rtracklayerTranscriptsFromEnsemblGtf(object)
    }



## FIXME RETHINK THIS APPROACH.
.standardizeGencodeToEnsembl <-
    function(object) {
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
.makeGenesFromRefSeqGff <-
    function(object) {
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
.makeGenesFromRefSeqGtf <-
    function(object) {
        assert(
            is(object, "GRanges"),
            isSubset(
                x = c("gene_biotype", "gene_id", "type"),
                y = colnames(mcols(object))
            ),
            areDisjointSets("gene_name", colnames(mcols(object)))
        )
        ## Define `gene_name` from `gene_id`.
        mcols(object)[["gene_name"]] <- mcols(object)[["gene_id"]]
        object <- object[mcols(object)[["type"]] == "gene"]
        object
    }



## Updated 2020-01-20.
.makeTranscriptsFromRefSeqGff <-
    function(object) {
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
.makeTranscriptsFromRefSeqGtf <-
    function(object) {
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



## FIXME RETHINK THIS APPROACH.
## This step ensures that `gene_id` and `transcript_id` columns are defined.
## Updated 2020-01-20.
.standardizeRefSeqToEnsembl <-
    function(object) {
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



## UCSC ========================================================================

## FIXME NEED TO ADD SUPPORT FOR THIS.

## .makeGenesFromUcscGtf
## .makeTranscriptsFromUcscGtf



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



.makeGenesFromWormBaseGTF <-
    function(object) {
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



.makeTranscriptsFromWormBaseGtf <-
    function(object) {
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
