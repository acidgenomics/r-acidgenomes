## RefSeq
##
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
.makeGenesFromRefSeqGFF3 <- function(object) {
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
.makeTranscriptsFromRefSeqGFF3 <- function(object) {
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
    cli_alert_info(sprintf(
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
