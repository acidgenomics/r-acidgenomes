## Updated 2021-01-27.
.makeGRangesFromRtracklayer <- function(
    file,
    level = c("genes", "transcripts"),
    ignoreVersion = FALSE,
    synonyms = FALSE
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
            object = paste("rtracklayer", provider, level, format),
            strict = TRUE
        )
    )
    tryCatch(
        expr = {
            what <- .getFun(funName)
        },
        error = function(e) {
            stop(sprintf("Unsupported GFF: '%s'.", basename(file)))
        }
    )
    gr <- do.call(what = what, args = list("object" = gr))
    metadata(gr) <- meta
    seqinfo <- .getSeqinfo(meta)
    if (is(seqinfo, "Seqinfo")) {
        seqinfo(gr) <- seqinfo[seqlevels(gr)]
    }
    .makeGRanges(
        object = gr,
        ignoreVersion = ignoreVersion,
        synonyms = synonyms
    )
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



## Updated 2021-01-27.
.rtracklayerEnsemblGenesGtf <-
    function(object) {
        assert(
            is(object, "GRanges"),
            isSubset(
                x = c(
                    "gene_id",
                    "gene_version",
                    "type"
                ),
                y = names(mcols(object))
            ),
            areDisjointSets(
                x = "gene_id_version",
                y = names(mcols(object))
            )
        )
        keep <- mcols(object)[["type"]] == "gene"
        object <- object[keep]
        assert(hasNoDuplicates(mcols(object)[["gene_id"]]))
        ## Match ensembldb versioned identifier convention.
        mcols(object)[["gene_id_version"]] <-
            paste(
                mcols(object)[["gene_id"]],
                mcols(object)[["gene_version"]],
                sep = "."
            )
        mcols(object)[["gene_version"]] <- NULL
        object
    }



## Updated 2021-01-27.
.rtracklayerEnsemblTranscriptsGtf <-
    function(object) {
        assert(
            is(object, "GRanges"),
            isSubset(
                x = c(
                    "gene_id",
                    "gene_version",
                    "transcript_id",
                    "transcript_version",
                    "type"
                ),
                y = names(mcols(object))
            ),
            areDisjointSets(
                x = c(
                    "gene_id_version",
                    "transcript_id_version"
                ),
                y = names(mcols(object))
            )
        )
        keep <- mcols(object)[["type"]] == "transcript"
        object <- object[keep]
        ## Match ensembldb versioned identifier convention.
        mcols(object)[["gene_id_version"]] <-
            paste(
                mcols(object)[["gene_id"]],
                mcols(object)[["gene_version"]],
                sep = "."
            )
        mcols(object)[["transcript_id_version"]] <-
            paste(
                mcols(object)[["transcript_id"]],
                mcols(object)[["transcript_version"]],
                sep = "."
            )
        mcols(object)[["gene_version"]] <- NULL
        mcols(object)[["transcript_version"]] <- NULL
        object
    }



## Updated 2021-01-27.
.rtracklayerEnsemblGenesGff <-
    function(object) {
        assert(
            is(object, "GRanges"),
            isSubset(
                x = c("Name", "biotype", "gene_id", "version"),
                y = names(mcols(object))
            ),
            areDisjointSets(
                x = c("gene_biotype", "gene_id_version", "gene_name"),
                y = names(mcols(object))
            )
        )
        keep <- !is.na(mcols(object)[["gene_id"]])
        object <- object[keep]
        assert(hasNoDuplicates(mcols(object)[["gene_id"]]))
        names(mcols(object))[
            names(mcols(object)) == "Name"] <- "gene_name"
        names(mcols(object))[
            names(mcols(object)) == "biotype"] <- "gene_biotype"
        ## Match ensembldb versioned identifier convention.
        mcols(object)[["gene_id_version"]] <-
            paste(
                mcols(object)[["gene_id"]],
                mcols(object)[["version"]],
                sep = "."
            )
        mcols(object)[["version"]] <- NULL
        object
    }



## Updated 2021-01-27.
.rtracklayerEnsemblTranscriptsGff <-
    function(object) {
        genes <- .rtracklayerEnsemblGenesGff(object)
        mcols(genes) <- removeNA(mcols(genes))
        assert(
            is(object, "GRanges"),
            isSubset(
                x = c("Name", "Parent", "biotype"),
                y = names(mcols(object))
            ),
            areDisjointSets(
                x = c("transcript_biotype", "transcript_name"),
                y = names(mcols(object))
            )
        )
        keep <- !is.na(mcols(object)[["transcript_id"]])
        object <- object[keep]
        assert(
            hasNoDuplicates(mcols(object)[["transcript_id"]]),
            allAreMatchingRegex(
                pattern = "^gene:",
                x = as.character(mcols(object)[["Parent"]])
            )
        )
        mcols(object) <- removeNA(mcols(object))
        names(mcols(object))[
            names(mcols(object)) == "biotype"] <- "transcript_biotype"
        names(mcols(object))[
            names(mcols(object)) == "Name"] <- "transcript_name"
        mcols(object)[["gene_id"]] <- gsub(
            pattern = "^gene:",
            replacement = "",
            x = as.character(mcols(object)[["Parent"]])
        )
        ## Match ensembldb versioned identifier convention.
        mcols(object)[["transcript_id_version"]] <-
            paste(
                mcols(object)[["transcript_id"]],
                mcols(object)[["version"]],
                sep = "."
            )
        mcols(object)[["version"]] <- NULL
        geneCols <- c(
            "gene_id",
            setdiff(
                x = names(mcols(genes)),
                y = names(mcols(object))
            )
        )
        mcols(object) <- leftJoin(
            x = mcols(object),
            y = mcols(genes)[geneCols],
            by = "gene_id"
        )
        object
    }



## FlyBase =====================================================================
## GTF:
## > [1] "source"            "type"              "score"
## > [4] "phase"             "gene_id"           "gene_symbol"
## > [7] "transcript_id"     "transcript_symbol" "#"



## Updated 2021-01-27.
.rtracklayerFlyBaseGenesGtf <-
    function(object) {
        assert(
            is(object, "GRanges"),
            isSubset(
                x = c("gene_id", "type"),
                y = names(mcols(object))
            )
        )
        keep <- mcols(object)[["type"]] == "gene"
        object <- object[keep]
        assert(hasNoDuplicates(mcols(object)[["gene_id"]]))
        object
    }



## Updated 2021-01-27.
.rtracklayerFlyBaseTranscriptsGtf <-
    function(object) {
        assert(
            is(object, "GRanges"),
            isSubset(
                x = c("transcript_id", "type"),
                y = names(mcols(object))
            )
        )
        keep <- grepl(
            pattern = paste(c("^pseudogene$", "RNA$"), collapse = "|"),
            x = mcols(object)[["type"]],
            ignore.case = TRUE
        )
        object <- object[keep]
        assert(hasNoDuplicates(mcols(object)[["transcript_id"]]))
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



## Updated 2021-01-27.
.rtracklayerGencodeGenesGtf <-
    function(object) {
        assert(
            is(object, "GRanges"),
            isSubset(
                x = c("gene_id", "type"),
                y = names(mcols(object))
            ),
            areDisjointSets(
                x = c("gene_id_version", "gene_version"),
                y = names(mcols(object))
            )
        )
        keep <- mcols(object)[["type"]] == "gene"
        object <- object[keep]
        assert(hasNoDuplicates(mcols(object)[["gene_id"]]))
        ## Match ensembldb versioned identifier convention.
        mcols(object)[["gene_id_version"]] <-
            mcols(object)[["gene_id"]]
        mcols(object)[["gene_id"]] <-
            stripGeneVersions(mcols(object)[["gene_id"]])
        object
    }



## Updated 2021-01-27.
.rtracklayerGencodeTranscriptsGtf <-
    function(object) {
        assert(
            is(object, "GRanges"),
            isSubset(
                x = c(
                    "gene_id",
                    "transcript_id",
                    "type"
                ),
                y = names(mcols(object))
            ),
            areDisjointSets(
                x = c(
                    "gene_id_version",
                    "gene_version",
                    "transcript_id_version",
                    "transcript_version"
                ),
                y = names(mcols(object))
            )
        )
        keep <- mcols(object)[["type"]] == "transcript"
        object <- object[keep]
        assert(hasNoDuplicates(mcols(object)[["transcript_id"]]))
        ## Match ensembldb versioned identifier convention.
        mcols(object)[["gene_id_version"]] <- mcols(object)[["gene_id"]]
        mcols(object)[["gene_id"]] <-
            stripGeneVersions(mcols(object)[["gene_id"]])
        mcols(object)[["transcript_id_version"]] <-
            mcols(object)[["transcript_id"]]
        mcols(object)[["transcript_id"]] <-
            stripGeneVersions(mcols(object)[["transcript_id"]])
        object
    }




## FIXME SIMPLIFY THE ASSERT HERE.
## FIXME WHAT ABOUT IGNORING VERSION HERE?

## Updated 2021-01-26.
.rtracklayerGenesFromGencodeGff <-
    function(object) {
        object <- .standardizeGencodeToEnsembl(object)
        assert(
            isSubset(
                x = c("ID", "gene_biotype", "gene_id", "gene_name", "type"),
                y = names(mcols(object))
            )
        )
        keep <- mcols(object)[["type"]] == "gene"
        object <- object[keep]
        ## Keep track of PAR genes (e.g. "ENSG00000002586.20_PAR_Y").
        mcols(object)[["gene_id"]] <- mcols(object)[["ID"]]
        mcols(object)[["ID"]] <- NULL
        object
    }



## Updated 2021-01-26.
.rtracklayerTranscriptsFromGencodeGff <-
    function(object) {
        ## FIXME SIMPLIFY THE ASSERT HERE.
        assert(
            isSubset(
                x = c(
                    "gene_biotype", "gene_id",
                    "transcript_biotype", "transcript_id"
                ),
                y = names(mcols(object))
            )
        )
        keep <- mcols(object)[["type"]] == "transcript"
        object <- object[keep]
        ## Keep track of PAR genes (e.g. "ENSG00000002586.20_PAR_Y").
        mcols(object)[["transcript_id"]] <- mcols(object)[["ID"]]
        mcols(object)[["ID"]] <- NULL
        object
    }



## RefSeq ======================================================================

## FIXME These parsers are currently experimental.
##       Using TxDb as current appraoch instead.

## FIXME NEED TO SANITIZE COLS FOR REFSEQ
## genes <- genes[!is.na(sanitizeNA(mcols(genes)[["gene_id"]]))]
## genes <- genes[is.na(sanitizeNA(mcols(genes)[["tx_id"]]))]



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



## Updated 2020-01-26.
.rtracklayerGenesFromRefSeqGff <-
    function(object) {
        assert(
            is(object, "GRanges"),
            isSubset(
                x = c(
                    "Dbxref",
                    "Name",
                    "gbkey",
                    "gene",
                    "gene_biotype",
                    "type"
                ),
                y = names(mcols(object))
            ),
            areDisjointSets(
                x = "gene_id",
                y = names(mcols(object))
            )
        )
        keep <- !is.na(mcols(object)[["gbkey"]])
        object <- object[keep]
        keep <- mcols(object)[["gbkey"]] == "Gene"
        object <- object[keep]


        ## Unpack the numeric gene identifier from the Dbxref column.
        ## FIXME RETHINK THIS APPROACH, DOESNT WORK WITH GTF...
        ## FIXME NEED TO MATCH CONVENTION USED FOR GTF, FOR SIMPLICITY.
        dbxref <- mcols(object)[["Dbxref"]]
        ids <- vapply(
            X = dbxref,
            FUN = function(x) {
                grep(pattern = "^GeneID:[0-9]+$", x = x, value = TRUE)
            },
            FUN.VALUE = character(1L)
        )

        ## FIXME THIS APROACH DOESNT WORK?
        ## NOT ALL ELEMENTS HAVE GENEID IN GTF FILE...

        ids <- gsub(pattern = "^GeneID:", replacement = "", x = ids)
        ids <- as.integer(ids)
        mcols(object)[["gene_id"]] <- ids
        names(mcols(object))[
            names(mcols(object)) == "gene"] <- "gene_name"
        ## FIXME NEED TO RENAME AND KEEP DBXREF KEY.
        object
    }



## Updated 2020-01-26.
.rtracklayerGenesFromRefSeqGtf <-
    function(object) {
        assert(
            is(object, "GRanges"),
            isSubset(
                x = c("db_xref", "gbkey", "gene_biotype", "gene_id", "type"),
                y = names(mcols(object))
            ),
            areDisjointSets("gene_name", names(mcols(object)))
        )
        keep <- mcols(object)[["gbkey"]] == "Gene"
        object <- object[keep]
        mcols(object)[["gene_id"]] <- mcols(object)[["gene"]]
        mcols(object)[["gene"]] <- NULL
        mcols(object)[["gene_name"]] <- mcols(object)[["gene_id"]]
        names(mcols(object))[
            names(mcols(object)) == "db_xref"] <- "dbxref"
        object
    }



## FIXME NEED TO ADD SUPPORT.

## Updated 2020-01-20.
.makeTranscriptsFromRefSeqGff <-
    function(object) {
        assert(
            is(object, "GRanges"),
            isSubset("Name", names(mcols(object))),
            areDisjointSets("transcript_name", names(mcols(object)))
        )
        ## Only keep annotations that map to `Name` column.
        keep <- !is.na(mcols(object)[["Name"]])
        object <- object[keep]
        object
    }



## FIXME NEED TO ADD SUPPORT.

## Updated 2020-01-20.
.makeTranscriptsFromRefSeqGtf <-
    function(object) {
        assert(
            is(object, "GRanges"),
            isSubset(c("transcript_id", "type"), names(mcols(object))),
            areDisjointSets(
                x = c("gene_name", "transcript_biotype", "transcript_name"),
                y = names(mcols(object))
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



## FIXME RETHINK THIS APPROACH, CONSIDER USING ABOVE.
## This step ensures that `gene_id` and `transcript_id` columns are defined.
## Updated 2020-01-20.
.standardizeRefSeqToEnsembl <-
    function(object) {
        assert(is(object, "GRanges"))
        mcols <- mcols(object)
        assert(
            isSubset(
                x = c("gene", "transcript_id"),
                y = names(mcols)
            ),
            areDisjointSets(
                x = c("gene_name", "transcript_name"),
                y = names(mcols)
            )
        )
        ## Ensure `gene_id` is defined.
        if (isTRUE(all(c("gene", "gene_id") %in% names(mcols)))) {
            ## Pick `gene_id` over `gene` column, if both are defined.
            ## This applies to GTF spec.
            keep <- setdiff(names(mcols), "gene")
            mcols <- mcols[keep]
        } else if ("gene" %in% names(mcols)) {
            ## Rename `gene` column to `gene_id`, matching Ensembl spec.
            ## This applies to GFF3 spec.
            names(mcols) <- sub("^gene$", "gene_id", names(mcols))
        }
        assert(isSubset(c("gene_id", "transcript_id"), names(mcols)))
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



## FIXME NEED TO CHECK THAT THIS WORKS.
## Updated 2021-01-15.
.rtracklayerGenesFromWormBaseGTF <-
    function(object) {
        assert(
            is(object, "GRanges"),
            isSubset(x = "gene_id", y = names(mcols(object))),
            areDisjointSets(x = "gene_name", y = names(mcols(object)))
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



## FIXME NEED TO CHECK THAT THIS WORKS.
## Updated 2021-01-15.
.rtracklayersTranscriptsFromWormBaseGtf <-
    function(object) {
        assert(
            is(object, "GRanges"),
            isSubset(
                x = c("gene_id", "transcript_id"),
                y = names(mcols(object))
            ),
            areDisjointSets(
                x = c("gene_name", "transcript_name"),
                y = names(mcols(object))
            )
        )
        ## Sanitize the `gene_id` column, which can contain some malformed
        ## entries prefixed with "Gene:". We want to keep these entries.
        mcols(object)[["gene_id"]] <- gsub(
            pattern = "^Gene:",
            replacement = "",
            x = mcols(object)[["gene_id"]]
        )
        ## Process using Ensembl conventions.
        object <- .makeTranscriptsFromEnsemblGTF(object)
        object
    }
