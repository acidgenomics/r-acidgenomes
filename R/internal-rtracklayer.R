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
        mcols(object)[["gene_id_version"]] <- mcols(object)[["gene_id"]]
        mcols(object)[["gene_id"]] <-
            stripGeneVersions(mcols(object)[["gene_id"]])
        mcols(object)[["transcript_id_version"]] <-
            mcols(object)[["transcript_id"]]
        mcols(object)[["transcript_id"]] <-
            stripGeneVersions(mcols(object)[["transcript_id"]])
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



## WormBase ====================================================================
##
## GTF:
## >  [1] "source"             "type"               "score"
## >  [4] "phase"              "gene_id"            "gene_source"
## >  [7] "gene_biotype"       "transcript_id"      "transcript_source"
## > [10] "transcript_biotype" "exon_number"        "exon_id"
## > [13] "protein_id"



## Updated 2021-01-27.
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



## Updated 2021-01-27.
.rtracklayersWormBaseTranscriptsGtf <-
    function(object) {
        assert(
            is(object, "GRanges"),
            isSubset(
                x = c("gene_id", "transcript_id", "type"),
                y = names(mcols(object))
            )
        )
        keep <- mcols(object)[["type"]] == "transcript"
        object <- object[keep]
        assert(
            hasNoDuplicates(mcols(object)[["transcript_id"]]),
            allAreMatchingRegex(
                pattern = "^WBGene[[:digit:]]{8}$",
                x = mcols(object)[["gene_id"]]
            )
        )
        object
    }
