## Updated 2021-08-05.
.makeGRangesFromRtracklayer <- function(
    file,
    level,
    ignoreVersion,
    synonyms
) {
    assert(isString(file))
    level <- match.arg(
        arg = level,
        choices = c("genes", "transcripts")
    )
    meta <- getGFFMetadata(file)
    assert(
        isString(meta[["format"]]),
        isString(meta[["provider"]])
    )
    meta[["level"]] <- level
    gr <- import(file = .cacheIt(file))
    assert(is(gr, "GRanges"))
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
        assert(
            any(keep),
            msg = "Failed to extract any genes."
        )
        object <- object[keep]
        assert(hasNoDuplicates(mcols(object)[["gene_id"]]))
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
        assert(
            any(keep),
            msg = "Failed to extract any transcripts."
        )
        object <- object[keep]
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
        assert(
            any(keep),
            msg = "Failed to extract any genes."
        )
        object <- object[keep]
        assert(hasNoDuplicates(mcols(object)[["gene_id"]]))
        names(mcols(object))[
            names(mcols(object)) == "Name"] <- "gene_name"
        names(mcols(object))[
            names(mcols(object)) == "biotype"] <- "gene_biotype"
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
        assert(
            any(keep),
            msg = "Failed to extract any transcripts."
        )
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
        assert(
            any(keep),
            msg = "Failed to extract any genes."
        )
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
        assert(
            any(keep),
            msg = "Failed to extract any transcripts."
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
        assert(
            any(keep),
            msg = "Failed to extract any genes."
        )
        object <- object[keep]
        assert(hasNoDuplicates(mcols(object)[["gene_id"]]))
        mcols(object)[["gene_id_version"]] <-
            mcols(object)[["gene_id"]]
        mcols(object)[["gene_id"]] <-
            stripGeneVersions(mcols(object)[["gene_id"]])
        object
    }



## Updated 2021-02-01.
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
        assert(
            any(keep),
            msg = "Failed to extract any transcripts."
        )
        object <- object[keep]
        mcols(object)[["transcript_id_version"]] <-
            mcols(object)[["transcript_id"]]
        mcols(object)[["transcript_id"]] <-
            stripTranscriptVersions(mcols(object)[["transcript_id"]])
        assert(hasNoDuplicates(mcols(object)[["transcript_id"]]))
        mcols(object)[["gene_id_version"]] <- mcols(object)[["gene_id"]]
        mcols(object)[["gene_id"]] <-
            stripGeneVersions(mcols(object)[["gene_id"]])
        object
    }



## Updated 2021-01-27.
.rtracklayerGencodeGenesGff <-
    function(object) {
        assert(
            is(object, "GRanges"),
            isSubset(
                x = c("ID", "gene_id", "type"),
                y = names(mcols(object))
            ),
            areDisjointSets(
                x = c(
                    "gene_id_no_version",
                    "gene_version"
                ),
                y = names(mcols(object))
            )
        )
        keep <- mcols(object)[["type"]] == "gene"
        assert(
            any(keep),
            msg = "Failed to extract any genes."
        )
        object <- object[keep]
        mcols(object)[["gene_id_version"]] <- mcols(object)[["ID"]]
        mcols(object)[["ID"]] <- NULL
        mcols(object)[["gene_id"]] <-
            stripGeneVersions(mcols(object)[["gene_id_version"]])
        assert(hasNoDuplicates(mcols(object)[["gene_id"]]))
        object
    }



## Updated 2021-01-27.
.rtracklayerGencodeTranscriptsGff <-
    function(object) {
        assert(
            isSubset(
                x = c(
                    "ID",
                    "gene_id",
                    "transcript_id",
                    "type"
                ),
                y = names(mcols(object))
            ),
            areDisjointSets(
                x = c(
                    "transcript_id_no_version",
                    "transcript_version"
                ),
                y = names(mcols(object))
            )
        )
        keep <- mcols(object)[["type"]] == "transcript"
        assert(
            any(keep),
            msg = "Failed to extract any transcripts."
        )
        object <- object[keep]
        mcols(object)[["transcript_id_version"]] <- mcols(object)[["ID"]]
        mcols(object)[["ID"]] <- NULL
        mcols(object)[["transcript_id"]] <-
            stripTranscriptVersions(mcols(object)[["transcript_id_version"]])
        assert(hasNoDuplicates(mcols(object)[["transcript_id"]]))
        mcols(object)[["gene_id_version"]] <-
            as.character(mcols(object)[["Parent"]])
        mcols(object)[["gene_id"]] <-
            stripGeneVersions(mcols(object)[["gene_id_version"]])
        object
    }



## RefSeq ======================================================================
## GTF:
##  [1] "source"             "type"               "score"
##  [4] "phase"              "gene_id"            "transcript_id"
##  [7] "db_xref"            "description"        "gbkey"
## [10] "gene"               "gene_biotype"       "pseudo"
## [13] "product"            "transcript_biotype" "exon_number"
## [16] "gene_synonym"       "model_evidence"     "tag"
## [19] "protein_id"         "note"               "exception"
## [22] "inference"          "anticodon"          "partial"
## [25] "transl_except"      "standard_name"      "codons"
## [28] "transl_table"

## GFF:
##  [1] "source"                    "type"
##  [3] "score"                     "phase"
##  [5] "ID"                        "Dbxref"
##  [7] "Name"                      "chromosome"
##  [9] "gbkey"                     "genome"
## [11] "mol_type"                  "description"
## [13] "gene"                      "gene_biotype"
## [15] "pseudo"                    "Parent"
## [17] "product"                   "transcript_id"
## [19] "gene_synonym"              "model_evidence"
## [21] "tag"                       "protein_id"
## [23] "Note"                      "experiment"
## [25] "function"                  "regulatory_class"
## [27] "standard_name"             "recombination_class"
## [29] "feat_class"                "rpt_type"
## [31] "rpt_unit_seq"              "exception"
## [33] "inference"                 "anticodon"
## [35] "partial"                   "start_range"
## [37] "end_range"                 "mobile_element_type"
## [39] "rpt_family"                "transl_except"
## [41] "satellite"                 "bound_moiety"
## [43] "Target"                    "assembly_bases_aln"
## [45] "assembly_bases_seq"        "bit_score"
## [47] "blast_aligner"             "blast_score"
## [49] "common_component"          "e_value"
## [51] "filter_score"              "for_remapping"
## [53] "gap_count"                 "hsp_percent_coverage"
## [55] "matchable_bases"           "matched_bases"
## [57] "num_ident"                 "num_mismatch"
## [59] "pct_coverage"              "pct_coverage_hiqual"
## [61] "pct_identity_gap"          "pct_identity_gapopen_only"
## [63] "pct_identity_ungap"        "rank"
## [65] "weighted_identity"         "lxr_locAcc_currStat_120"
## [67] "not_for_annotation"        "Gap"
## [69] "consensus_splices"         "exon_identity"
## [71] "identity"                  "idty"
## [73] "matches"                   "product_coverage"
## [75] "splices"                   "map"
## [77] "part"                      "merge_aligner"
## [79] "lxr_locAcc_currStat_35"    "direction"
## [81] "rpt_unit_range"            "exon_number"
## [83] "number"                    "allele"
## [85] "align_id"                  "batch_id"
## [87] "crc32"                     "curated_alignment"
## [89] "promoted_rank"             "qtaxid"
## [91] "Is_circular"               "country"
## [93] "isolation-source"          "note"
## [95] "tissue-type"               "codons"
## [97] "transl_table"



## Updated 2021-08-06.
.rtracklayerRefSeqGenesGtf <-
    function(object) {
        assert(
            is(object, "GRanges"),
            isSubset(
                x = c(
                    "gene",
                    "gene_id",
                    "type"
                ),
                y = names(mcols(object))
            )
        )
        keep <- mcols(object)[["type"]] == "gene"
        assert(
            any(keep),
            msg = "Failed to extract any genes."
        )
        object <- object[keep]
        names(mcols(object))[
            names(mcols(object)) == "gene_id"] <- "parent_gene_id"
        names(mcols(object))[
            names(mcols(object)) == "gene"] <- "gene_id"
        assert(hasNoDuplicates(mcols(object)[["parent_gene_id"]]))
        object
    }



## Updated 2021-08-06.
.rtracklayerRefSeqTranscriptsGtf <-
    function(object) {
        genes <- .rtracklayerRefSeqGenesGtf(object)
        genesMcols <- mcols(genes)[
            ,
            c(
                "parent_gene_id",
                "gene_biotype",
                "description"
            ),
            drop = FALSE
        ]
        keep <- complete.cases(genesMcols)
        genesMcols <- genesMcols[keep, , drop = FALSE]
        assert(
            hasNoDuplicates(genesMcols[["parent_gene_id"]]),
            is(object, "GRanges"),
            isSubset(
                x = c(
                    "transcript_id",
                    "type"
                ),
                y = names(mcols(object))
            )
        )
        keep <- mcols(object)[["type"]] == "transcript"
        assert(
            any(keep),
            msg = "Failed to extract any transcripts."
        )
        object <- object[keep]
        assert(hasNoDuplicates(mcols(object)[["transcript_id"]]))
        names(mcols(object))[
            names(mcols(object)) == "gene_id"] <- "parent_gene_id"
        names(mcols(object))[
            names(mcols(object)) == "gene"] <- "gene_id"
        cols <- c(
            setdiff(
                x = colnames(mcols(object)),
                y = colnames(genesMcols)
            ),
            "parent_gene_id"
        )
        mcols(object) <- mcols(object)[, cols]
        mcols(object) <- leftJoin(
            x = mcols(object),
            y = genesMcols,
            by = "parent_gene_id"
        )
        object
    }



## Updated 2021-08-06.
.rtracklayerRefSeqGenesGff <-
    function(object) {
        assert(
            is(object, "GRanges"),
            isSubset(
                x = c("ID", "Parent", "description", "gene"),
                y = names(mcols(object))
            ),
            areDisjointSets(
                x = "gene_id",
                y = names(mcols(object))
            )
        )
        keep <-
            !is.na(mcols(object)[["gbkey"]]) &
            mcols(object)[["gbkey"]] == "Gene"
        assert(
            any(keep),
            msg = "Failed to extract any genes."
        )
        object <- object[keep]
        names(mcols(object))[names(mcols(object)) == "ID"] <- "parent_gene_id"
        assert(hasNoDuplicates(mcols(object)[["parent_gene_id"]]))
        names(mcols(object))[ names(mcols(object)) == "gene"] <- "gene_id"
        assert(all(grepl(
            pattern = "^gene-",
            x = mcols(object)[["parent_gene_id"]]
        )))
        mcols(object)[["parent_gene_id"]] <- gsub(
            pattern = "^gene-",
            replacement = "",
            x = mcols(object)[["parent_gene_id"]]
        )
        object
    }



## Updated 2021-08-06.
.rtracklayerRefSeqTranscriptsGff <-
    function(object) {
        genes <- .rtracklayerRefSeqGenesGff(object)
        genesMcols <- mcols(genes)[
            ,
            c("parent_gene_id", "gene_biotype", "description"),
            drop = FALSE
        ]
        keep <- complete.cases(genesMcols)
        genesMcols <- genesMcols[keep, , drop = FALSE]
        assert(
            hasNoDuplicates(genesMcols[["parent_gene_id"]]),
            is(object, "GRanges"),
            isSubset(
                x = c("Parent", "gene", "transcript_id"),
                y = names(mcols(object))
            ),
            areDisjointSets(
                x = "gene_id",
                y = names(mcols(object))
            )
        )
        keep <- !is.na(mcols(object)[["transcript_id"]])
        assert(
            any(keep),
            msg = "Failed to extract any transcripts."
        )
        object <- object[keep]
        ## Only keep transcript annotations that map to a parent gene.
        ## This is somewhat slow and may be optimizable.
        keep <- bapply(
            X = mcols(object)[["Parent"]],
            FUN = function(x) {
                any(grepl(pattern = "^gene-", x = x))
            }
        )
        assert(
            any(keep),
            msg = "Failed to match transcripts against parent genes."
        )
        object <- object[keep]
        ## e.g. "NM_000218.3".
        assert(allAreMatchingRegex(
            pattern = "^[A-Z]{2}_[0-9]+\\.[0-9]+$",
            x = mcols(object)[["transcript_id"]]
        ))
        ## Ensure that matching transcripts contain a unique gene parent.
        assert(
            all(bapply(
                X = mcols(object)[["Parent"]],
                FUN = isScalar
            )),
            msg = "Elements do not contain a unique gene parent."
        )
        mcols(object)[["parent_gene_id"]] <- vapply(
            X = mcols(object)[["Parent"]],
            FUN = function(x) {
                sub(pattern = "^gene-", replacement = "", x = x[[1L]])
            },
            FUN.VALUE = character(1L),
            USE.NAMES = FALSE
        )
        mcols(object)[["Parent"]] <- NULL
        names(mcols(object))[names(mcols(object)) == "gene"] <- "gene_id"
        cols <- c(
            setdiff(
                x = colnames(mcols(object)),
                y = colnames(genesMcols)
            ),
            "parent_gene_id"
        )
        mcols(object) <- mcols(object)[, cols]
        mcols(object) <- leftJoin(
            x = mcols(object),
            y = genesMcols,
            by = "parent_gene_id"
        )
        object
    }



## UCSC ========================================================================
## GTF:
## [1] "source"        "type"          "score"         "phase"
## [5] "gene_id"       "transcript_id" "gene_name"     "exon_number"
## [9] "exon_id"
##
## These files are complicated to parse, and we're simply handing off to TxDb
## in the current version of the package.



## WormBase ====================================================================
## GTF (2021-08-05):
##  [1] "source"             "type"               "score"
##  [4] "phase"              "gene_id"            "gene_version"
##  [7] "gene_source"        "gene_biotype"       "gene_name"
## [10] "transcript_id"      "transcript_source"  "transcript_biotype"
## [13] "exon_number"        "exon_id"            "protein_id"



## Updated 2021-08-05.
.rtracklayerWormBaseGenesGtf <-
    function(object) {
        assert(
            is(object, "GRanges"),
            isSubset(
                x = c("gene_id", "type"),
                y = names(mcols(object))
            )
        )
        keep <- mcols(object)[["type"]] == "gene"
        assert(
            any(keep),
            msg = "Failed to extract any transcripts."
        )
        object <- object[keep]
        assert(
            hasNoDuplicates(mcols(object)[["gene_id"]]),
            allAreMatchingRegex(
                pattern = "^WBGene[[:digit:]]{8}$",
                x = mcols(object)[["gene_id"]]
            )
        )
        object
    }



## Updated 2021-08-05.
.rtracklayerWormBaseTranscriptsGtf <-
    function(object) {
        assert(
            is(object, "GRanges"),
            isSubset(
                x = c("gene_id", "transcript_id", "type"),
                y = names(mcols(object))
            )
        )
        ## Keep track of gene-level metadata, that we'll join below.
        genes <- .rtracklayerWormBaseGenesGtf(object)
        geneCols <- grep(
            pattern = "^gene_",
            x = colnames(mcols(genes)),
            value = TRUE
        )
        geneMcols <- mcols(genes, use.names = FALSE)[, geneCols]
        keep <- mcols(object)[["type"]] == "transcript"
        assert(
            any(keep),
            msg = "Failed to extract any transcripts."
        )
        object <- object[keep]
        assert(
            hasNoDuplicates(mcols(object)[["transcript_id"]]),
            allAreMatchingRegex(
                pattern = "^WBGene[[:digit:]]{8}$",
                x = mcols(object)[["gene_id"]]
            )
        )
        mcols <- mcols(object)
        cols <- c(
            setdiff(
                x = colnames(mcols),
                y = colnames(geneMcols)
            ),
            "gene_id"
        )
        mcols <- mcols[, cols]
        mcols <- leftJoin(x = mcols, y = geneMcols, by = "gene_id")
        mcols(object) <- mcols
        object
    }
