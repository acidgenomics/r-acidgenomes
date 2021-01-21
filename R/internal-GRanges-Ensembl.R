## Ensembl
##
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
.makeGenesFromEnsemblGFF3 <- function(object) {
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



.makeTranscriptsFromEnsemblGFF3 <- function(object) {
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
