# Release notes

## AcidGenomes 0.6.2 (2023-11-21)

Major changes:

- `EnsemblToNcbi` and `NcbiToEnsembl`: Added support for disabling `strict`
  mode, which is useful for mapping all genes in a reference genome. This is
  not a breaking change, as strict mode remains enabled by default. Also
  reworked the internal code to speed up 1:1 mapping return.

Minor changes:

- `Hgnc`: Renamed `"symbol"` column to `"geneName"` and `"name"` to
  `"description"`, better matching naming conventions in other functions.

## AcidGenomes 0.6.1 (2023-10-13)

Minor changes:

- Renamed `HumanToMouse` to `JaxHumanMouse`, which indicates the source
  (Jackson Laboratory) more clearly.
- Bug fixes for GENCODE and RefSeq genome downloads.

## AcidGenomes 0.6.0 (2023-10-04)

New functions / classes:

- `HumanToMouse`: Downloads human-to-mouse gene mappings from the Jackson
  Laboratory MGI server. Returns with unique 1:1 mappings by default, which
  can be disabled with the `unique` argument.

Major changes:

- Renamed some identifier classes to use `"To"` instead of the numeric `"2"` in
  the class name. Sorry Prince, but this just looks weird for some functions.
  Applies to `Ensembl2Ncbi`, `Gene2Symbol`, `Ncbi2Ensembl`, `Protein2Gene`,
  and `Tx2Gene`.
- Hardened class checks via `validObject` for our custom classes, ensuring
  that metadata is consistently slotted, including date and package version.
- Switched to HTTPS downloads internally for some functions, as the Ensembl
  FTP server has been prone to timeouts lately. This does not affect our
  primary genome downloading functions, which still use FTP for speed.

Minor changes:

- Migrated internal data from `inst/extdata` to save internally inside the
  package via `sysdata.rda` file. Applies to `detectOrganism` and
  `mapNcbiTaxId` currently.

## AcidGenomes 0.5.1 (2023-08-01)

Major changes:

- `downloadGencodeGenome`: This genome download function now sanitizes the
  transcriptome FASTA to only include transcript identifiers in the header
  without additional information separated by the `"|"` (pipe) delimeter. This
  approach is not commonly used in FASTA files, and results in unwanted
  downstream behavior when quantifying at transcript level using kallisto and
  minimap2. Note that salmon can currently handle this edge case when setting
  the `--gencode` flag during the genome index step. This action is
  non-destructive and returns a "transcriptome_fixed" FASTA file. We are now
  symlinking this fixed file by default, but the unmodified original is
  retained in the transcriptome download folder.

Minor changes:

- `currentEnsemblGenomeBuild`: Fixed internal REST API query to Ensembl server.
  This now requires `"content-type=application/json"` to be defined in the URL,
  otherwise the Ensembl server returns text instead of JSON.
- `currentEnsemblVersion`: Now parses `current_README` file on FTP server
  instead of the top level `README`. We changed this because the `README`
  symlink can break during Ensembl release updates (e.g. 109 to 110, in
  progress).

## AcidGenomes 0.5.0 (2023-04-27)

Starting a new release series to denote potential breaking changes with legacy
objects saved with `entrezId` instead of `ncbiGeneId`.

New functions:

- Added utility functions to easily map gene symbols to stable gene identifiers
  at Ensembl, HGNC, and NCBI (Entrez): `mapGeneNamesToEnsembl`,
  `mapGeneNamesToHGNC`, `mapGeneNamesToNCBI`. These are covered against
  _Homo sapiens_ and _Mus musculus_.

Major changes:

- Classes that extend `GRanges` (e.g. `EnsemblGenes`, `EnsemblTranscripts`,
  `GencodeGenes`, `GencodeTranscripts`) now intentionally fail class checks if
  `entrezId` is defined instead of `ncbiGeneId` in `mcols` metadata. This makes
  downstream handoff to GSEA functions in AcidGSEA easier to manage. For
  legacy objects, use `updateObject` to resolve this check.
- `makeGRangesFromEnsembl` and `makeGRangesFromGFF` now attempt to fetch
  additional useful gene metadata, including gene synonyms from the Ensembl FTP
  server when applicable. This currently applies to gene annotation files from
  Ensembl and GENCODE. Note that extra metadata is not supported for legacy
  _Homo sapiens_ GRCh37 genome build.

Minor changes:

- Renamed `EntrezGeneInfo` function to `NcbiGeneInfo`.
- Removed `geneSynonyms` from `NAMESPACE`. Consider using `NcbiGeneInfo` or
  `HGNC` instead for synonym information.
- Removed `HGNC2Ensembl` and `MGI2Ensembl`. Just use `HGNC` and `MGI` function
  return instead.
- `downloadEnsemblGenome` now downloads additional useful metadata files.
- Added `updateObject` support to update legacy objects that may fail new
  `entrezId` class checks.

## AcidGenomes 0.4.8 (2023-02-16)

Minor changes:

- `makeTx2GeneFromGFF`, `Tx2Gene`: Ensure that rownames are defined for RefSeq
  genome annotations, which are constructed from `GenomicRangesList` method.

## AcidGenomes 0.4.7 (2023-02-10)

New functions:

- `gencodeReleaseHistory`: This function scrapes the GENCODE website to return
  the full release history for either human or mouse genomes.

Minor changes:

- `currentEnsemblVersion`: Fix for breaking change on Ensembl FTP server.
  The file this function parses has been renamed from `current_README` to
  simply `README`.
- `mapGencodeToEnsembl`: Now using `gencodeReleaseHistory` internally to
  dynamically fetch metadata directly from the GENCODE website, instead of
  relying on an internal CSV mapping file. This helps avoid having to update
  the package every time a new Ensembl/GENCODE release comes out.

## AcidGenomes 0.4.6 (2023-02-09)

Minor changes:

- Migrated `requireNamespaces` import from AcidBase to goalie.
- Updated dependency versions.
- Now allowing Ensembl `Seqinfo` fetch step to fail, to avoid unit test
  issues with _Mus musulus_ Ensembl 90 GTF file (for bcbioRNASeq).

## AcidGenomes 0.4.5 (2023-01-31)

Minor changes:

- `downloadUCSCGenome`: Added a manual override for defaulting to `"hg38"`
  for _Homo sapiens_, which has switched over to experimental `"hs1"` T2T
  genome build.
- `mapGencodeToEnsembl`: Added support for mapping _Mus musculus_ releases.
- Internal `Seqinfo` generators now default to pulling annotations from Ensembl
  for GENCODE reference, rather than UCSC. Note that UCSC `Seqinfo` function
  is currently broken for hg38, but will be fixed in pending GenomeInfoDb
  v1.34.8.
- Updated unit tests to check against latest genome build releases.

## AcidGenomes 0.4.4 (2022-10-25)

Minor changes:

- Compatibility update to provide support for breaking changes introduced with
  pipette 0.10.0 release.
- `EntrezGeneInfo`: Added column name checks. Renamed `xTaxId` to `taxonomyId`.
- Migrated all `importFrom` calls to `imports.R` file.
- Removed BiocParallel bplapply usage.
- `export`: Updated methods to match new conventions defined in pipette.
- Improved support for quick detection of NCBI taxonomic group for commonly
  used model organisms. Previously human and mouse were supported. Added support
  for zebrafish, rat, worm, fruitfly, and yeast.

## AcidGenomes 0.4.3 (2022-06-09)

Minor chagnes:

- `makeGRangesFromEnsembl`: Hardened internal code to suppress spurious warnings
  from rtracklayer due to masking of `download.file` function. See
  [issue #71](https://github.com/lawremi/rtracklayer/issues/71) for details.
- Improved temporary directory creation and deletion using `tempdir2` and
  `unlink2` internally, which improves support for continuous integration (CI)
  checks on Windows.

## AcidGenomes 0.4.2 (2022-05-27)

Minor changes:

- `Ensembl2Entrez`: Simplified `format` `"1:1"` handling, based on new code
  approach used in AcidGSEA package. This now keeps track of original rownames,
  which is necessary for bcbioRNASeq clusterProfiler R Markdown template.

## AcidGenomes 0.4.1 (2022-05-24)

Minor changes:

- Updated lintr checks and testthat unit tests.

## AcidGenomes 0.4.0 (2022-05-04)

Major changes:

- S4 classes containing identifier mappins now are set to contain `DFrame`
  instead of `DataFrame`, as the previous definition approach no longer works
  with Bioconductor 3.15. This applies to `Ensembl2Entrez`, `Entrez2Ensembl`,
  `EntrezGeneInfo`, `Gene2Symbol`, `HGNC`, `HGNC2Ensembl`, `Protein2Gene`,
  and `Tx2Gene` classes.
- `downloadGencodeGenome`: Added support for Entrez and RefSeq identifiers,
  which are now defined in `mcols` of `GRanges` objects.
- Reworked internal code of rtracklayer runner functions. In particular,
  improved identifier handling for RefSeq GFF and GTF runners, adding support
  for extraction of Entrez gene identifiers.

Minor changes:

- Reformatted code using styler conventions.
- Removed internal dependencies on stringr, in favor of base R when possible
  and stringi package otherwise.
- Hardened import of external files using readr package, where applicable.
  Applies primarily to Entrez, HGNC, and MGI file import.
- `makeTx2GeneFromFASTA`: Added `ignoreVersions` argument, which is now enabled
  by default, to match other `Tx2Gene` functions.
- Updated unit tests to cover most recent versions of genome assembly files.

## AcidGenomes 0.3.0 (2022-03-11)

New functions:

- `mapGencodeToEnsembl`: Convenience function for mapping human GENCODE release
  (e.g. `39`) to corresponding Ensembl release (e.g. `105`).

Minor changes:

- Reworked NAMESPACE to inherit from Bioconductor packages, migrating from
  previous approach that used generics defined only in AcidGenerics.
- Relaxed `GenomicRanges` `mcols` column name validity checks to no longer
  require strict camelCase formatting. This check can run into issues when
  slotting into a `DESeqDataSet` object, which appends non-camelCase-formatted
  columns into the `mcols` (corresponding to the `rowRanges` of the object).
- Migrated some generics to AcidGenerics, and reexporting here (non-breaking
  change): `Ensembl2Entrez`, `Entrez2Ensembl`, `Gene2Symbol`, `Tx2Gene`.
- Updated references from `GRanges` to `GenomicRanges` (virtual class)
  where applicable in the package.
- Split out documentation on `Ensembl2Entrez` and `Entrez2Ensembl` methods
  into separate files.
- Increased the verbosity of `setMethod` calls, particular in the `signature`
  argument, where applicable.
- Calls to `AnnotationHub` internally now should never prompt the user about
  whether to create cache directory. This is achieved by setting `ask = FALSE`
  internally.
- Updated internal class checks inside of `assert` calls to check for
  `GenomicRanges` virtual class rather than `GRanges`.
- Improved alphabetical sorting of exported S4 methods.
- Reexporting some additional functions from GenomeInfoDb: `Seqinfo`,
  `genome`, `seqinfo`, and `seqlevels`. Note that corresponding assignment
  methods (if applicable) are intentionally not reexported here.
- Genome download functions no longer cache using BiocFileCache by default.

## AcidGenomes 0.2.20 (2022-01-08)

Minor changes:

- `Tx2Gene` class check: disabling check that looks for identical transcript
  and gene identifiers. This check is not compatible with the
  _Saccharomyces cerevisiae_ (sacCer3) reference genome. Thanks for pointing
  this out @amizeranschi.

## AcidGenomes 0.2.19 (2021-09-13)

Minor changes:

- `getEnsDb` / `makeGRangesFromEnsembl`: Quieted down package loading from
  Bioconductor when obtaining annotations for GRCh37 (`EnsDb.Hsapiens.v75`
  release package).
- Improved CLI messages to use S4 class when applicable.

## AcidGenomes 0.2.18 (2021-09-03)

Minor changes:

- Improved CLI message formatting.
- `downloadRefSeqGenome`, `downloadUCSCGenome`: Improved the `genomeBuild`
  documentation, with more specific examples.

## AcidGenomes 0.2.17 (2021-08-11)

Major changes:

- `Gene2Symbol`: Hardened internal identifier mapping code in `switch` call
  to support `format` argument. Improved unit testing for expected behavior
  of `format` argument. Fixed "1:1" mapping to split based on `geneName` column
  rather than `geneId` column.
- `Tx2Gene`: Improved code coverge and cleaned up internal `complete.cases`
  handling.
- All identifier mapping classes (e.g. `Ensembl2Entrez`, `Gene2Symbol`,
  `Tx2Gene`) now check for `complete.cases` in S4 validity methods.

Minor changes:

- `currentEnsemblVersion` and `mapHumanOrthologs` working examples are now
  re-enabled, wrapped in a `try` call.

## AcidGenomes 0.2.16 (2021-08-09)

Major changes:

- Removed `mapEnsemblBuildToUCSC` and `mapUCSCBuildToEnsembl` functions. Also
  removed mapping support for UCSC genome build names (e.g. "hg38") inside of
  `makeGRangesFromEnsembl` calls, since this is not technically the correct
  genome build name.
- `downloadEnsemblGenome`, `downloadGencodeGenome`, etc. now support file
  caching by default with `cache = TRUE` argument.

Minor changes:

- Improved genome download code coverage.
- Split out `stripGeneVersions` and `stripTranscriptVersions` documentation
  into separate files.
- Reorganized S4 method documentation to be alphabetical consistently.
- `export`: Hardened `Tx2Gene` method to ensure that rownames are consistently
  removed prior to export. Noticed that this was an issue with UCSC genome
  build download.

## AcidGenomes 0.2.15 (2021-07-27)

Minor changes:

- Added support for new _Mus musculus_ GRCm39 genome build. Updated internal
  code for `mapUCSCBuildToEnsembl` and `downloadEnsemblGenome`, in particular.
  Note that `*_chr_patch_hapl_scaff` GFF and GTF files are no longer available
  on the Egnyte FTP server for GRCm39 (only GRCm38 and GRCh38).

## AcidGenomes 0.2.14 (2021-06-10)

Minor changes:

- `Gene2Symbol`: Improve handling when gene identifiers are integer, such as
  is the case with NCBI Entrez gene identifiers.

## AcidGenomes 0.2.13 (2021-05-18)

- Internal fixes to provide compatibility for R 4.1 release.
- Updated dependencies to support Bioconductor 3.13.

## AcidGenomes 0.2.12 (2021-04-27)

Major changes:

- Removed some Bioconductor packages from imports: AnnotationDbi, AnnotationHub,
  GenomeInfoDb, and ensembldb. This helps reduce package loading time and avoid
  unwanted BiocManager messages from appearing at startup (due to loading of
  AnnotationHub). These are included as suggested packages, which should not
  be problematic, since they are frequently used.

## AcidGenomes 0.2.11 (2021-03-19)

Minor changes:

- `HGNC` now returns columns with split values as `CharacterList`, instead of
  as character strings containing "|".

## AcidGenomes 0.2.10 (2021-03-15)

Minor changes:

- `mapHumanOrthologs`: Hardened mouse-to-human matching.
- `makeGRangesFromEnsembl`: No longer hard-coding minimum release version check
  at 87, in case older releases are ported to AnnotationHub in a future release.
- Revert back to `ignoreVersion = TRUE` by default for genome annotation
  importers, as this is typically what users expect by default.

## AcidGenomes 0.2.9 (2021-03-03)

Minor changes:

- `Gene2Symbol` functions now preserve metadata, as expected. This was
  causing pointillism package to error, due to unwanted breaking change.
- `Tx2Gene`: Improved consistency of metadata return, ensuring `call` and
  `synonyms` are not defined.
- Renamed internal "acidGenomes" metadata key to "packageVersion", for
  consistency with conventions used in other Acid Genomics packages.

## AcidGenomes 0.2.8 (2021-03-02)

Minor changes:

- Relaxed validity checks for `EnsemblGenes` and `EnsemblTranscripts`.

## AcidGenomes 0.2.7 (2021-02-26)

Minor changes:

- `makeGRangesFromGFF`: Improved support and code coverage for handling of
  bcbio-nextgen `ref-transcripts.gtf` genome file.

## AcidGenomes 0.2.6 (2021-02-25)

Minor changes:

- Updated basejump dependency versions.
- `EntrezGeneInfo`: Improved column formatting.

## AcidGenomes 0.2.5 (2021-02-17)

Minor changes:

- Genome downloader functions (e.g. `downloadEnsemblGenome`) now return
  relative symlinks instead of absolute paths.
- Bug fix for `mapHumanOrthologs` internal join step. Now returns `humanGeneId`
  and `humanGeneName` columns instead of `hgncId` and `hgncName` columns, which
  technically were incorrect, since these map to Ensembl.

## AcidGenomes 0.2.4 (2021-02-13)

New functions:

- `EntrezGeneInfo`: New utility for obtaining gene annotations from NCBI.

Major changes:

- `geneSynonyms`: Reworked internal code, extending `EntrezGeneInfo`.

## AcidGenomes 0.2.3 (2021-02-10)

Minor changes:

- Reverted back to using vroom as importer for `HGNC` and `MGI2Ensembl`.

## AcidGenomes 0.2.2 (2021-02-09)

Minor changes:

- Reduced the number of reexported functions.

## AcidGenomes 0.2.1 (2021-02-08)

Minor changes:

- Now including some reexports from GenomicRanges and IRanges.

## AcidGenomes 0.2.0 (2021-02-02)

New functions:

- Added new genome download functions, migrated from previous Python approach
  defined in py-koopa package.
- Switched from cli to AcidCLI for interactive messages.

Major changes:

- In metadata columns, renamed "ID" to "Id", for stricter lower camel case
  formatting.
- Now ensuring GenomicRanges gets attached as a dependency.
- Now using BiocFileCache (via `pipette::cacheURL` internally) to automatically
  cache GFF/GTF files when used in `makeGRangesFromGFF`.
- Renamed `ignoreTxVersion` to simply `ignoreVersion`, where applicable.
  We want this setting to also apply at gene level.
- GFF, TxDb, and ensembldb parser functions now default to
  `ignoreVersion = FALSE`. Previous releases of AcidGenomes and basejump had
  this set to `ignoreVersion = TRUE` by default. Note that both modes are now
  non-destructive.
- `GRanges` `mcols` now return with `tx` prefix instead of `transcript`.
- `GRanges` `mcols` now use strict camel case formatting
  (e.g. `geneId` instead of `geneID`).

## AcidGenomes 0.1.1 (2020-10-12)

Minor changes:

- Relaxed stringency of internal organism name checks. Applies to
  `makeGene2SymbolFromEnsembl` for example, which is causing running examples
  to fail in pointillism without a fix.
- Updated minimum dependency versions.

## AcidGenomes 0.1.0 (2020-10-07)

Initial release, consisting of functions migrated from basejump.
