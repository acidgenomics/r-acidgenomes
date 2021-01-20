#' Get the metadata (directives) from a GFF file
#'
#' @export
#' @note Updated 2021-01-20.
#'
#' @inheritParams AcidRoxygen::params
#'
#' @details
#' Matches lines beginning with `#!<key> <value>` or `##<key>: <value>`
#'
#' @section GFF3:
#' Lines beginning with '##' are directives (sometimes called pragmas or
#' meta-data) and provide meta-information about the document as a whole. Blank
#' lines should be ignored by parsers and lines beginning with a single '#' are
#' used for human-readable comments and can be ignored by parsers. End-of-line
#' comments (comments preceded by # at the end of and on the same line as a
#' feature or directive line) are not allowed.
#'
#' @seealso
#' - https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
#'
#' @return `DataFrame` or `NULL`.
#'
#' @examples
#' ## > getGFFMetadata("Homo_sapiens.GRCh38.102.gtf.gz")
#' ## DataFrame with 5 rows and 2 columns
#' ##                      key                 value
#' ##              <character>           <character>
#' ## 1           genome-build            GRCh38.p13
#' ## 2         genome-version                GRCh38
#' ## 3            genome-date               2013-12
#' ## 4 genome-build-accession NCBI:GCA_000001405.28
#' ## 5 genebuild-last-updated               2020-09
#'
#' ## > getGFFMetadata("Homo_sapiens.GRCh38.102.gff3.gz")
#' ## DataFrame with 200 rows and 2 columns
#' ##                        key                 value
#' ##                <character>           <character>
#' ## 1              gff-version                     3
#' ## 2          sequence-region         1 1 248956422
#' ## 3          sequence-region        10 1 133797422
#' ## 4          sequence-region        11 1 135086622
#' ## 5          sequence-region        12 1 133275309
#' ## ...                    ...                   ...
#' ## 196           genome-build            GRCh38.p13
#' ## 197         genome-version                GRCh38
#' ## 198            genome-date               2013-12
#' ## 199 genome-build-accession NCBI:GCA_000001405.28
#' ## 200 genebuild-last-updated               2020-09
#'
#' ## > getGFFMetadata("Homo_sapiens.GRCh37.87.gtf.gz")
#' ## DataFrame with 5 rows and 2 columns
#' ##                      key                 value
#' ##              <character>           <character>
#' ## 1           genome-build            GRCh37.p13
#' ## 2         genome-version                GRCh37
#' ## 3            genome-date               2009-02
#' ## 4 genome-build-accession NCBI:GCA_000001405.14
#' ## 5 genebuild-last-updated               2013-09
#'
#' ## > getGFFMetadata("Homo_sapiens.GRCh37.87.gff3.gz")
#' ## DataFrame with 90 rows and 2 columns
#' ##                        key                 value
#' ##                <character>           <character>
#' ## 1              gff-version                     3
#' ## 2          sequence-region         1 1 249250621
#' ## 3          sequence-region        10 1 135534747
#' ## 4          sequence-region        11 1 135006516
#' ## 5          sequence-region        12 1 133851895
#' ## ...                    ...                   ...
#' ## 86            genome-build            GRCh37.p13
#' ## 87          genome-version                GRCh37
#' ## 88             genome-date               2009-02
#' ## 89  genome-build-accession NCBI:GCA_000001405.14
#' ## 90  genebuild-last-updated               2013-09
#'
#' ## > getGFFMetadata("gencode.v36.annotation.gtf.gz")
#' ## DataFrame with 5 rows and 2 columns
#' ##           key                  value
#' ##   <character>            <character>
#' ## 1 description evidence-based annot..
#' ## 2    provider                GENCODE
#' ## 3     contact gencode-help@ebi.ac.uk
#' ## 4      format                    gtf
#' ## 5        date             2020-09-30
#'
#' ## > getGFFMetadata("gencode.v36.annotation.gff3.gz")
#' ## DataFrame with 31 rows and 2 columns
#' ##                 key                  value
#' ##         <character>            <character>
#' ## 1       gff-version                      3
#' ## 2       description evidence-based annot..
#' ## 3          provider                GENCODE
#' ## 4           contact gencode-help@ebi.ac.uk
#' ## 5            format                   gff3
#' ## ...             ...                    ...
#' ## 27  sequence-region       chr21 1 46709983
#' ## 28  sequence-region       chr22 1 50818468
#' ## 29  sequence-region       chrX 1 156040895
#' ## 30  sequence-region        chrY 1 57227415
#' ## 31  sequence-region           chrM 1 16569
#'
#' ## > getGFFMetadata("gencode.v36lift37.annotation.gtf.gz")
#' ## DataFrame with 5 rows and 2 columns
#' ##           key                  value
#' ##   <character>            <character>
#' ## 1 description evidence-based annot..
#' ## 2    provider                GENCODE
#' ## 3     contact gencode-help@ebi.ac.uk
#' ## 4      format                   gff3
#' ## 5        date             2020-09-30
#'
#' ## > getGFFMetadata("gencode.v36lift37.annotation.gff3.gz")
#' ## DataFrame with 31 rows and 2 columns
#' ##                 key                  value
#' ##         <character>            <character>
#' ## 1       gff-version                      3
#' ## 2       description evidence-based annot..
#' ## 3          provider                GENCODE
#' ## 4           contact gencode-help@ebi.ac.uk
#' ## 5            format                   gff3
#' ## ...             ...                    ...
#' ## 27  sequence-region       chr21 1 48129895
#' ## 28  sequence-region       chr22 1 51304566
#' ## 29  sequence-region       chrX 1 155270560
#' ## 30  sequence-region        chrY 1 59373566
#' ## 31  sequence-region           chrM 1 16569
#'
#' ## > getGFFMetadata("gencode.vM25.annotation.gtf.gz")
#' ## DataFrame with 5 rows and 2 columns
#' ##           key                  value
#' ##   <character>            <character>
#' ## 1 description evidence-based annot..
#' ## 2    provider                GENCODE
#' ## 3     contact gencode-help@ebi.ac.uk
#' ## 4      format                    gtf
#' ## 5        date             2020-03-24
#'
#' ## > getGFFMetadata("gencode.vM25.annotation.gff3.gz")
#' ## DataFrame with 28 rows and 2 columns
#' ##                 key                  value
#' ##         <character>            <character>
#' ## 1       gff-version                      3
#' ## 2       description evidence-based annot..
#' ## 3          provider                GENCODE
#' ## 4           contact gencode-help@ebi.ac.uk
#' ## 5            format                   gff3
#' ## ...             ...                    ...
#' ## 24  sequence-region       chr18 1 90702639
#' ## 25  sequence-region       chr19 1 61431566
#' ## 26  sequence-region       chrX 1 171031299
#' ## 27  sequence-region        chrY 1 91744698
#' ## 28  sequence-region           chrM 1 16299
#'
#' ## > getGFFMetadata("GCF_000001405.38_GRCh38.p12_genomic.gtf.gz")
#' ## DataFrame with 4 rows and 2 columns
#' ##                      key                  value
#' ##              <character>            <character>
#' ## 1            gtf-version                    2.2
#' ## 2           genome-build             GRCh38.p12
#' ## 3 genome-build-accession NCBI_Assembly:GCF_00..
#' ## 4      annotation-source NCBI Homo sapiens An..
#'
#' ## > getGFFMetadata("GCF_000001405.38_GRCh38.p12_genomic.gff.gz")
#' ## DataFrame with 1194 rows and 2 columns
#' ##                         key                  value
#' ##                 <character>            <character>
#' ## 1               gff-version                      3
#' ## 2          gff-spec-version                   1.21
#' ## 3                 processor       NCBI annotwriter
#' ## 4              genome-build             GRCh38.p12
#' ## 5    genome-build-accession NCBI_Assembly:GCF_00..
#' ## ...                     ...                    ...
#' ## 1190                species https://www.ncbi.nlm..
#' ## 1191        sequence-region   NT_113949.2 1 177381
#' ## 1192                species https://www.ncbi.nlm..
#' ## 1193        sequence-region    NC_012920.1 1 16569
#' ## 1194                species https://www.ncbi.nlm..
#'
#' ## > getGFFMetadata("dmel-all-r6.37.gtf.gz")
#' ## NULL
#'
#' ## > getGFFMetadata("dmel-all-r6.37.gff.gz")
#' ## DataFrame with 1874 rows and 2 columns
#' ##                   key                  value
#' ##           <character>            <character>
#' ## 1         gff-version                      3
#' ## 2             species http://www.ncbi.nlm...
#' ## 3    feature-ontology ftp://ftp.flybase.or..
#' ## 4        genome-build          FlyBase r6.37
#' ## 5     sequence-region 211000022279464 1 1412
#' ## ...               ...                    ...
#' ## 1870  sequence-region 211000022278407 1 1191
#' ## 1871  sequence-region 211000022280133 1 2162
#' ## 1872  sequence-region 211000022279437 1 1433
#' ## 1873  sequence-region 211000022280009 1 1425
#' ## 1874  sequence-region 211000022278174 1 1154
#'
#' ## > getGFFMetadata("c_elegans.PRJNA13758.WS279.canonical_geneset.gtf.gz")
#' ## DataFrame with 1 row and 2 columns
#' ##                 key       value
#' ##         <character> <character>
#' ## 1 genebuild-version       WS279
#'
#' ## > getGFFMetadata("c_elegans.PRJNA13758.WS279.annotations.gff3.gz")
#' ## DataFrame with 20 rows and 2 columns
#' ##                 key          value
#' ##         <character>    <character>
#' ## 1       gff-version              3
#' ## 2   sequence-region   I 1 15072434
#' ## 3   sequence-region  II 1 15279421
#' ## 4   sequence-region III 1 13783801
#' ## 5   sequence-region  IV 1 17493829
#' ## ...             ...            ...
#' ## 16  sequence-region   I 1 15072434
#' ## 17  sequence-region  II 1 15279421
#' ## 18  sequence-region  IV 1 17493829
#' ## 19  sequence-region   X 1 17718942
#' ## 20  sequence-region   V 1 20924180
getGFFMetadata <- function(file, nMax = Inf) {
    file <- .cacheIt(file)
    lines <- import(
        file = file,
        format = "lines",
        nMax = nMax,
        quiet = TRUE
    )
    pattern <- "^(#!|#+)([a-z-]+)(:)?\\s(.+)$"
    lines <- grep(pattern = pattern, x = lines, value = TRUE)
    if (!hasLength(lines)) return(NULL)
    mat <- str_match(
        string = grep(pattern = pattern, x = lines, value = TRUE),
        pattern = pattern
    )
    assert(is.matrix(mat), hasRows(mat))
    df <- as(mat, "DataFrame")
    df <- df[, c(3L, 5L), drop = FALSE]
    colnames(df) <- c("key", "value")
    df
}
